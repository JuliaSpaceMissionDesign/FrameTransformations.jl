export format_camelcase, format_snakecase

function format_camelcase(::Type{T}, s::S) where {T,S<:AbstractString}
    words = split(s, r"[\_, \s, \-]")
    if length(words) > 1
        T(join(uppercasefirst.(lowercase.(words))))
    elseif length(words) == 1
        T(uppercasefirst(s))
    else
        throw(error("$s cannot be formatted in CamelCase."))
    end
end

function format_snakecase(::Type{T}, s::S) where {T,S<:AbstractString}
    return T(join(lowercase.(split(replace(s, r"[\-\.\s]" => "_"), r"(?=[A-Z])")), "_"))
end
