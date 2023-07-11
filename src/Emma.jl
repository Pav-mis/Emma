using Serialization
using Artifacts
using BioSequences
using FASTX
using XGBoost
using Logging
using Unicode

const emmamodels = "Emma/emma_vertebrate_models"

include("circularity.jl")
include("feature.jl")
include("orfs.jl")
include("tRNAs.jl")
include("rRNAs.jl")
include("gff.jl")
include("visuals.jl")

const minORF = 150

function remove_starts_in_tRNAs!(starts, codons, tRNAs, glength)
    fms = getproperty.(tRNAs, :fm)
    for (startvector, codonvector) in zip(starts, codons)
        startsintRNA = []
        for (i,s) in enumerate(startvector)
            if any(circularin.(s, fms, glength))
                push!(startsintRNA, i)
            end
        end
        deleteat!(startvector, startsintRNA)
    end
end

function add_stops_at_tRNAs!(stops, tRNAs, glength)
    for t in tRNAs
        push!(stops[mod1(t.fm.target_from, 3)], t.fm.target_from)
        minus1 = mod1(t.fm.target_from-1, glength)
        push!(stops[mod1(minus1, 3)], minus1)
        minus2 = mod1(t.fm.target_from-2, glength)
        push!(stops[mod1(minus2, 3)], minus2)
    end
    sort!(stops[1])
    sort!(stops[2])
    sort!(stops[3])
    return stops
end

function main(infile::String, outfile::String; svgfile="", loglevel="Info")

    global_logger(ConsoleLogger(loglevel == "debug" ? Logging.Debug : Logging.Info))

    @info "$infile"
    target = FASTA.Record()
    reader = open(FASTA.Reader, infile)
    read!(reader, target)
    id = FASTA.identifier(target)
    genome = CircularSequence(FASTA.sequence(LongDNA{4}, target))
    glength = length(genome)
    rev_genome = reverse_complement(genome)
    @info "$id\t$glength bp"

    #extend genome
    extended_genome = genome[1:glength+100]
    writer = open(FASTA.Writer, "tmp.extended.fa")
    write(writer, FASTA.Record(id, extended_genome))
    close(writer)

    #find tRNAs
    trn_matches = parse_trn_alignments(cmsearch("trn", "all_trn.cm"), glength)
    @debug trn_matches
    filter!(x -> x.fm.evalue < 1e-5, trn_matches)
    filter!(x -> x.fm.target_from <= glength, trn_matches)
    ftrns = get_best_trns(filter(x->x.fm.strand=='+',trn_matches), glength)
    rtrns = get_best_trns(filter(x->x.fm.strand=='-',trn_matches), glength)
    #check for overlapping tRNAs
    overlapped = get_overlapped_trns(sort(ftrns, by=x->x.fm.target_from), glength)
    append!(overlapped, get_overlapped_trns(sort(rtrns, by=x->x.fm.target_from), glength))
    trn_matches = append!(ftrns, rtrns)
    #for overlapped tRNAs, generate polyadenylated version
    for (cma, trunc_end) in overlapped
        trnseq = cma.fm.strand == '+' ? genome.sequence[cma.fm.target_from:trunc_end] : rev_genome.sequence[cma.fm.target_from:trunc_end]
        trnseq_polyA = trnseq * dna"AAAAAAAAAA"
        polyA_matches = parse_trn_alignments(cmsearch(cma.fm.query, "trn", trnseq_polyA), 0)
        isempty(polyA_matches) && continue
        trn_match = polyA_matches[1]
        if trn_match.fm.evalue < cma.fm.evalue
            #construct modified CMA with new evalue, new end point, polyA length
            newcma = tRNA(FeatureMatch(cma.fm.id, cma.fm.query, cma.fm.strand, trn_match.fm.model_from, trn_match.fm.model_to, cma.fm.target_from, 
                circulardistance(cma.fm.target_from, trunc_end, glength) + 1, trn_match.fm.evalue),
                trn_match.anticodon, (trn_match.fm.target_length) - (trunc_end - cma.fm.target_from))
            #delete old match from trn_matches
            deleteat!(trn_matches, findfirst(x -> x == cma, trn_matches))
            #add new CMA
            push!(trn_matches, newcma)
        end
    end
    @debug trn_matches
    @info "found $(length(trn_matches)) tRNA genes"

    #find rRNAs
    #search for rrns using hmmsearch
    rrns = parse_tbl(rrnsearch(), glength)
    @debug rrns
    #fix ends using flanking trn genes
    ftRNAs = filter(x->x.fm.strand == '+', trn_matches)
    sort!(ftRNAs, by=x->x.fm.target_from)
    rtRNAs = filter(x->x.fm.strand == '-', trn_matches)
    sort!(rtRNAs, by=x->x.fm.target_from)
    fix_rrn_ends!(rrns, ftRNAs, rtRNAs, glength)
    @debug rrns

    #find CDSs
    fstarts, fstartcodons = getcodons(genome, startcodon)
    remove_starts_in_tRNAs!(fstarts, fstartcodons, ftRNAs, glength)
    fstops = codonmatches(genome, stopcodon)
    add_stops_at_tRNAs!(fstops, ftRNAs, glength)
    @debug fstops

    rstarts, rstartcodons = getcodons(rev_genome, startcodon)
    remove_starts_in_tRNAs!(rstarts, rstartcodons, rtRNAs, glength)
    rstops = codonmatches(rev_genome, stopcodon)
    add_stops_at_tRNAs!(rstops, rtRNAs, glength)

    cds_matches = parse_domt(orfsearch(id, genome, fstarts, fstops, rstarts, rstops, minORF), glength)    
    @debug cds_matches
    #fix start & stop codons
    #load XGBoost model
    startcodon_model = Booster(DMatrix[],model_file=joinpath(emmamodels, "xgb.model"))
    fhmms = filter(x->x.strand == '+', cds_matches)
    rhmms = filter(x->x.strand == '-', cds_matches)
    fix_start_and_stop_codons!(fhmms, ftRNAs, fstarts, fstartcodons, fstops, startcodon_model, glength)
    fix_start_and_stop_codons!(rhmms, rtRNAs, rstarts, rstartcodons, rstops, startcodon_model, glength)

    @info "found $(length(fhmms) + length(rhmms)) protein-coding genes"

    gffs = writeGFF(outfile, id, genome, rev_genome, append!(fhmms,rhmms), trn_matches, rrns)
    #drawgenome(svgfile, id, glength, gffs)
end

if abspath(PROGRAM_FILE) == @__FILE__ #running as script
    if isdir(ARGS[1]) && isdir(ARGS[2])
        fafiles = filter!(x->endswith(x,".fa") || endswith(x,".fasta"), readdir(ARGS[1], join=true))
        for fasta in fafiles
            accession = first(split(basename(fasta),"."))
            main(fasta, joinpath(ARGS[2], accession * ".gff"))
        end
    elseif length(ARGS) == 2
        main(ARGS[1], ARGS[2])
    elseif length(ARGS) == 3
        if startswith(ARGS[3], "svgfile")
            main(ARGS[1], ARGS[2]; svgfile=last(split(ARGS[3], "=")))
        elseif startswith(ARGS[3], "loglevel")
            main(ARGS[1], ARGS[2]; loglevel=last(split(ARGS[3], "=")))
        end
    elseif length(ARGS) == 4
        svgfile::String
        loglevel::String
        if startswith(ARGS[3], "svgfile")
            svgfile=last(split(ARGS[3], "="))
        elseif startswith(ARGS[3], "loglevel")
            loglevel=last(split(ARGS[3], "="))
        end
        if startswith(ARGS[4], "svgfile")
            svgfile=last(split(ARGS[4], "="))
        elseif startswith(ARGS[4], "loglevel")
            loglevel=last(split(ARGS[4], "="))
        end
    end
end

#ARGS[1] = fasta input
#ARGS[2] = gff output
#ARGS[3] = svg output
