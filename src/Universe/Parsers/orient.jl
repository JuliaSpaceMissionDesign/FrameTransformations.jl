using Basic.Orient: parse_iauconstants, generate_iauangles!

function _parse_orient!(gen::String, bid::Int, name, data::D1) where {D1<:AbstractDict}

    # IAU orient constants
    iauconst = parse_iauconstants([bid], data[:constants])

    @info "[Universe/Orient] Generating code for $name ($bid)..."
    # generate IAU angles
    gen = generate_iauangles!(gen, bid, iauconst)
    return gen
end
