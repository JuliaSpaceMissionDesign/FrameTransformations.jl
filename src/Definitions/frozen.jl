"""
    add_axes_frozen!(fr, name, id, frozen, e)

Create a frozen axes version of `frozen` at epoch `e`, with `name` and `id` and add 
it to `frames`.

!!! warning 
    The parent axes are automatically assigned to the `frozen` parent.
"""
function add_axes_frozen!(
    fr::FrameSystem{O, N, S}, name::Symbol, id::Int, frozen, epoch::Epoch{S2}
) where {O, N, S, S2}

    eS = convert(S(), epoch)
    fixid = axes_id(fr, frozen) 

    if !(has_axes(fr, fixid))
        throw(
            ErrorException(
                "no axes with ID $frozen is registered in the given frame system"
            )
        )
    end
    
    mid = get_mappedid(axes_graph(fr), fixid)
    @show node = get_mappednode(axes_graph(fr), mid)
    parid = node.parentid

    # Compute rotation 
    dcm = rotation3(fr, parid, fixid, j2000s(eS))[1]

    # Create new axes
    return add_axes_fixedoffset!(fr, name, id, parid, dcm)
end