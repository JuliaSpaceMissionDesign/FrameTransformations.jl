export parse_iso

function tryparsenext_base10(
    str::AbstractString, i::Int, len::Int, min_width::Int=1, max_width::Int=0
)
    i > len && return nothing
    min_pos = min_width <= 0 ? i : i + min_width - 1
    max_pos = max_width <= 0 ? len : min(i + max_width - 1, len)
    d::Int64 = 0
    @inbounds while i <= max_pos
        c, ii = iterate(str, i)::Tuple{Char,Int}
        if '0' <= c <= '9'
            d = d * 10 + (c - '0')
        else
            break
        end
        i = ii
    end
    if i <= min_pos
        return nothing
    else
        return d, i
    end
end

"""
    parse_iso(s::S) where {S<: AbstractString}

Parse ISO datetime string.

### Input 

- `s` -- ISO datetime `String` in the format `YYYY-MM-DDThh:mm:ss.ffffffff`

### Output 

A `Tuple` containing `year`, `month`, `day`, `hour`, `minute`, `second` and 
`millisecond` is parsed if the string is ISO otherwise an error is throw.
"""
function parse_iso(s::S) where {S<:AbstractString}
    i, end_pos = firstindex(s), lastindex(s)

    local dy
    dm = dd = Int64(1)
    th = tm = ts = tms = Int64(0)

    let val = tryparsenext_base10(s, i, end_pos, 1)
        val === nothing && @goto error
        dy, i = val
        i > end_pos && @goto done
    end

    c, i = iterate(s, i)::Tuple{Char,Int}
    c != '-' && @goto error
    i > end_pos && @goto done

    let val = tryparsenext_base10(s, i, end_pos, 1, 2)
        val === nothing && @goto error
        dm, i = val
        i > end_pos && @goto done
    end

    c, i = iterate(s, i)::Tuple{Char,Int}
    c != '-' && @goto error
    i > end_pos && @goto done

    let val = tryparsenext_base10(s, i, end_pos, 1, 2)
        val === nothing && @goto error
        dd, i = val
        i > end_pos && @goto done
    end

    c, i = iterate(s, i)::Tuple{Char,Int}
    c != 'T' && @goto error
    i > end_pos && @goto done

    let val = tryparsenext_base10(s, i, end_pos, 1, 2)
        val === nothing && @goto error
        th, i = val
        i > end_pos && @goto done
    end

    c, i = iterate(s, i)::Tuple{Char,Int}
    c != ':' && @goto error
    i > end_pos && @goto done

    let val = tryparsenext_base10(s, i, end_pos, 1, 2)
        val === nothing && @goto error
        tm, i = val
        i > end_pos && @goto done
    end

    c, i = iterate(s, i)::Tuple{Char,Int}
    c != ':' && @goto error
    i > end_pos && @goto done

    let val = tryparsenext_base10(s, i, end_pos, 1, 2)
        val === nothing && @goto error
        ts, i = val
        i > end_pos && @goto done
    end

    c, i = iterate(s, i)::Tuple{Char,Int}
    c != '.' && @goto error
    i > end_pos && @goto done

    let val = tryparsenext_base10(s, i, end_pos, 1, 8)
        val === nothing && @goto error
        tms, j = val
        tms *= 10.0^(-(j - i))
    end

    @label done
    return (dy, dm, dd, th, tm, ts, Float64(tms))

    @label error
    throw(
        ArgumentError(
            "[Tempo] Invalid ISO Epoch string! " *
            "The correct format is YYYY-MM-DDThh:mm:ss.ffff!",
        ),
    )
end
