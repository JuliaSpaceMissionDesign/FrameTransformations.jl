
"""
    add_axes_frozen!(frames, name, id, frozen, e)

Create a frozen axes version of `frozen` at epoch `e`, with `name` and `id` and add 
it to `frames`.

!!! warning 
    The parent axes are automatically assigned to the `frozen` parent.
"""
function add_axes_frozen!(
    frames::FrameSystem{O, N, S}, name::Symbol, id::Int, frozen, epoch::Epoch{S2}
) where {O, N, S, S2}

    eS = convert(S(), epoch)
    fixid = axes_id(frames, frozen) 

    if !(has_axes(frames, fixid))
        throw(
            ErrorException(
                "no axes with ID $frozen is registered in the given frame system"
            )
        )
    end
    
    mid = get_mappedid(axes_graph(frames), fixid)
    @show node = get_mappednode(axes_graph(frames), mid)
    parid = node.parentid

    # Compute rotation 
    dcm = rotation3(frames, parid, fixid, j2000s(eS))[1]

    # Create new axes
    return add_axes_fixedoffset!(frames, name, id, parid, dcm)
end