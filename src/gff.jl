struct GFF
    source::String
    ftype::String
    fstart::String
    fend::String
    score::String
    strand::Char
    phase::String
    attributes::String
end

function FMcoords2GFF(strand::Char, start::Integer, length::Integer, glength::Integer)
    gffstart = strand == '+' ? start : mod1(reverse_complement(start + length - 1, glength), glength)
    gffend = strand == '+' ? start + length - 1 : reverse_complement(start, glength)
    if gffend < gffstart; gffend += glength; end
    gffstart, gffend
end

function FMcoords2GFF(fm::FeatureMatch, glength::Integer)
    FMcoords2GFF(fm.strand, fm.target_from, fm.target_length, glength::Integer)
end

function CDS2GFF(cds::FeatureMatch, genome::CircularSequence, rev_genome::CircularSequence, trns::Vector{tRNA})
    glength = length(genome)
    attributes = "Name=" * cds.query
    cdsstart = cds.target_from
    cdsstop = cds.target_from + cds.target_length - 1 + 3
    trnidx = findfirst(t->circularin(t.fm.target_from, cdsstop-2, 3, glength), trns)
    if !isnothing(trnidx)
        stopcodon = cds.strand == '+' ? copy(genome[cdsstop-2:cdsstop-2+circulardistance(cdsstop-2, trns[trnidx].fm.target_from, glength)-1]) : copy(rev_genome[cdsstop-2:cdsstop-2+circulardistance(cdsstop-2, trns[trnidx].fm.target_from, glength)-1])
        cdsstop = trns[trnidx].fm.target_from - 1
        while length(stopcodon) < 3
            push!(stopcodon, DNA_A)
        end
        attributes *= ";Note=putative $stopcodon stop codon is completed by the addition of 3' A residues to the mRNA"
    end
    gffstart, gffend = FMcoords2GFF(cds.strand, cdsstart, circulardistance(cdsstart, cdsstop+1, glength), glength)
    return GFF("Emma", "CDS", string(gffstart), string(gffend), string(cds.evalue), cds.strand, "0", attributes)
end

function tRNA2GFF(trn::tRNA, glength::Integer)
    gffstart, gffend = FMcoords2GFF(trn.fm, glength)
    attributes = "Name=" * trn.fm.query * "-" * trn.anticodon
    if trn.polyA > 0
        attributes *= ";Note=tRNA completed by post-transcriptional addition of " * string(trn.polyA)
        attributes *= trn.polyA > 1 ? " As" : " A"
    end
    if typeof(attributes) != Missing
        return GFF("Emma", "tRNA", string(gffstart), string(gffend), string(trn.fm.evalue), trn.fm.strand, ".", attributes)
    else 
        return nothing
    end
end

function rRNA2GFF(rrn::FeatureMatch, glength::Integer)
    gffstart, gffend = FMcoords2GFF(rrn, glength)
    attributes = "Name=" * rrn.query
    return GFF("Emma", "rRNA", string(gffstart), string(gffend), string(rrn.evalue), rrn.strand, ".", attributes)
end

function writeGFF(outfile::String, id::AbstractString, genome::CircularSequence, rev_genome::CircularSequence, cds_matches::Vector{FeatureMatch},
             trn_matches::Vector{tRNA}, rRNAs::Vector{FeatureMatch})

    gffs = GFF[]
    genome_length = length(genome)

    function writeone(out::IO, gff::Union{Nothing, GFF})
        if typeof(gff) != Nothing
            write(out, join([id, gff.source, gff.ftype,gff.fstart,gff.fend,gff.score,gff.strand,gff.phase,gff.attributes], "\t"), "\n")
            push!(gffs, gff)
        end
    end

    open(outfile, "w") do out
        for cds in cds_matches
            gff = CDS2GFF(cds, genome, rev_genome, trn_matches)
            writeone(out, gff)
        end
        for trn in trn_matches
            gff = tRNA2GFF(trn, genome_length)
            writeone(out, gff)
        end
        for rrn in rRNAs
            gff = rRNA2GFF(rrn, genome_length)
            writeone(out, gff)
        end
    end
    return gffs
end
struct GFF
    source::String
    ftype::String
    fstart::String
    fend::String
    score::String
    strand::Char
    phase::String
    attributes::String
end

function HMMmatch2GFF(cds::HMMmatch, genome_length::Integer)
    cdsstart = cds.ali_from
    cdsstop = cds.ali_to + 2
    attributes = "Name=" * cds.query
    stopphase = mod(cdsstop-cdsstart+1, 3)
    if stopphase ≠ 0
        As_to_add = 3-stopphase
        attributes *= ";Note=transcript completed by post-transcriptional addition of " * string(As_to_add)
        attributes *= As_to_add > 1 ? " As" : " A"
    end
    if cds.strand == '-'
        tmp = cdsstart
        cdsstart = reverse_complement(cdsstop, genome_length)
        cdsstop = reverse_complement(tmp, genome_length)
    end
    if cdsstart > genome_length
        cdsstart = mod1(cdsstart, genome_length)
        cdsstop = mod1(cdsstop, genome_length)
    end
    return GFF("Emma", "CDS", string(cdsstart), string(cdsstop), string(cds.Evalue), cds.strand, "0", attributes)
end

function CMAlignment2GFF(trn::CMAlignment_trn, glength::Integer)
    startstring = trn.tstrand =='+' ? string(trn.tfrom) : string(reverse_complement(trn.tto, glength))
    finishstring = trn.tstrand =='+' ? string(trn.tto) : string(reverse_complement(trn.tfrom, glength))
    attributes = "Name=" * trn.query * "-" * trn.anticodon
    if trn.polyA > 0
        attributes *= ";Note=tRNA completed by post-transcriptional addition of " * string(trn.polyA)
        attributes *= trn.polyA > 1 ? " As" : " A"
    end
    if typeof(attributes) != Missing
        return GFF("Emma", "tRNA", startstring, finishstring, string(trn.Evalue), trn.tstrand, ".", attributes)
    else 
        return nothing
    end
end

function rRNA2GFF(rrn::rRNA, glength::Integer)
    start, stop = gettermini(rrn, glength)
    attributes = "Name=" * rrn.stop[1].query
    evalue = getevalue(rrn)
    return GFF("Emma", "rRNA", string(start), string(stop), string(evalue), rrn.stop[1].tstrand, ".", attributes)
end

function writeGFF(outfile::String, id::String, genome_length::Integer, cds_matches::Vector{HMMmatch},
             trn_matches::Vector{CMAlignment_trn}, rRNAs::NamedTuple{(:rrnL, :rrnS), Tuple{rRNA, rRNA}})

    gffs = GFF[]

    function writeone(out::IO, gff::GFF)
        write(out, join([id, gff.source, gff.ftype,gff.fstart,gff.fend,gff.score,gff.strand,gff.phase,gff.attributes], "\t"), "\n")
        push!(gffs, gff)
    end

    open(outfile, "w") do out
        for cds in cds_matches
            gff = HMMmatch2GFF(cds, genome_length)
            writeone(out, gff)
        end
        for trn in trn_matches
            gff = CMAlignment2GFF(trn, genome_length)
            if typeof(gff) != Nothing
                writeone(out, gff)
            end
        end
        gff = rRNA2GFF(rRNAs.rrnL, genome_length)
        writeone(out, gff)
        gff = rRNA2GFF(rRNAs.rrnS, genome_length)
        writeone(out, gff)
    end

    return gffs
end
