
# unified struct to cover all matches of HMMs or CMs; target coordinates are from the 5' end of the strand.
struct FeatureMatch
    id::String
    query::String
    strand::Char
    model_from::Int
    model_to::Int
    target_from::Int
    target_length::Int
    evalue::Float64
end

function circularin(x::Integer, f::FeatureMatch, c::Integer)
    circularin(x, f.target_from, f.target_length, c)
end

function circularoverlap(f1::FeatureMatch, f2::FeatureMatch, c::Integer)
    if circularin(f1.target_from, f2, c)
        return min(f1.target_from + f1.target_length, f2.target_from + f2.target_length) - f1.target_from
    elseif circularin(f2.target_from, f1, c)
        return min(f1.target_from + f1.target_length, f2.target_from + f2.target_length) - f2.target_from
    end
    return 0
end

function circulardistance(m1::FeatureMatch, m2::FeatureMatch, c::Integer)
    return circulardistance(m1.target_from, m2.target_from, c)
end

function merge_matches(m1::FeatureMatch, m2::FeatureMatch, glength::Integer)
    return FeatureMatch(m1.id, m1.query, m1.strand, min(m1.model_from, m2.model_from), max(m1.model_to, m2.model_to), m1.target_from,
        max(m1.target_length, circulardistance(m1.target_from, m2.target_from + m2.target_length, glength)), min(m1.evalue, m2.evalue))
end

# merge adjacent matches to same model
function rationalise_matches!(matches::Vector{FeatureMatch}, glength::Integer)
    length(matches) < 2 && return matches
    for i in matches, j in matches
        i == j && continue
        i.strand ≠ j.strand && continue
        # if not same gene and substantially overlap, only keep the best one
        if i.query ≠ j.query
            if circularoverlap(i, j, glength) > 0.5 * min(i.target_length, j.target_length)
                @debug i, j, circularoverlap(i, j, glength), 0.5 * min(i.target_length, j.target_length)
                todelete = i.evalue > j.evalue ? i : j
                deleteat!(matches, findfirst(x->x==todelete,matches))
                return(rationalise_matches!(matches, glength))
            else
                continue
            end
        end
        modeldistance = j.model_from - i.model_to
        modellength = max(i.model_to, j.model_to) - min(i.model_from, j.model_from)
        matchdistance = circulardistance(i.target_from + i.target_length - 1, j.target_from, glength)
        tolerance = 0.1
        if circularin(j.target_from, i, glength) || (matchdistance - modeldistance)^2 < tolerance * modellength^2 # arbitrary tolerance
            merged_match = merge_matches(i, j, glength)
            deleteat!(matches, findfirst(x->x==i,matches))
            deleteat!(matches, findfirst(x->x==j,matches))
            push!(matches, merged_match)
            return(rationalise_matches!(matches, glength))
        end
    end
    matches
end

