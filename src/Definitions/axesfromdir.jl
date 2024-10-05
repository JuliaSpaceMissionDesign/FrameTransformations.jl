"""
    function add_axes_twodir!(frames::FrameSystem{O,T}, name::Symbol, id, parent, 
        dir1::Symbol, dir2::Symbol, seq::Symbol; project::Bool=false) where {O,T}

Add a set of axes to `frames` based on two directions. 

A right-handed coordinate system is generated based on the specified sequence direction (`seq`), 
which determines the order in which the vectors are used to define the basis. 
The `project` flag specifies whether the resulting axes are inertial or not.

### See also 
See also [`add_axes_projected!`](@ref) and [`add_axes_rotating!`](@ref).
"""
function add_axes_twodir!(
    frames::FrameSystem{O,T}, name::Symbol, id, parent, dir1::Symbol, dir2::Symbol, seq::Symbol;
    project::Bool=false
) where {O,T}

    # Check directions 
    if !(has_direction(frames, dir1))
        throw(
            ArgumentError("No direction with name $dir1 available.")
        )
    end

    if !(has_direction(frames, dir2))
        throw(
            ArgumentError("No direction with name $dir2 available.")
        )
    end

    if !(has_axes(frames, parent))
        throw(
            ArgumentError("No axes with id $pid available.")
        )
    end

    fun = t -> twodir_to_dcm(
        direction3(frames, dir1, parent, t), direction3(frames, dir2, parent, t), seq
    )

    if project
        add_axes_projected!(frames, name, id, parent, fun)
    else
        add_axes_rotating!(frames, name, id, parent, fun)
    end
end

function _two_dir_basis(a::AbstractVector, b::AbstractVector,
    seq::Symbol, fc::Function)

    if seq == :XY
        w = fc(a, b)
        v = fc(w, a)
        u = a

    elseif seq == :YX
        w = fc(b, a)
        u = fc(a, w)
        v = a

    elseif seq == :XZ
        v = fc(b, a)
        w = fc(a, v)
        u = a

    elseif seq == :ZX
        v = fc(a, b)
        u = fc(v, a)
        w = a

    elseif seq == :YZ
        u = fc(a, b)
        w = fc(u, a)
        v = a

    elseif seq == :ZY
        u = fc(b, a)
        v = fc(a, u)
        w = a
    else
        throw(ArgumentError("Invalid rotation sequence $seq."))
    end

    return u, v, w
end

function twodir_to_dcm(a::AbstractVector, b::AbstractVector, seq::Symbol)
    ut, vt, wt = _two_dir_basis(a, b, seq, cross3)
    u, v, w = unitvec(ut), unitvec(vt), unitvec(wt)

    @inbounds dcm = DCM((u[1], v[1], w[1], u[2], v[2], w[2], u[3], v[3], w[3]))
    return dcm
end