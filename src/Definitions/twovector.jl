
function add_axes_twovectors!(
    frames::FrameSystem{O, N}, name::Symbol, id::Int, parentid::Int,
    from1::Int, to1::Int, from2::Int, to2::Int, seq::Symbol; 
    inertial::Bool=false
) where {O, N}

    if from1 == to1 || from2 == to2
        throw(
            ArgumentError("from/to points shall have a different value")
        )
    end

    if inertial
        funs = FrameAxesFunctions{O, N}(
            t -> Rotation{O}(
                twovectors_to_dcm(
                    vector3(frames, from1, to1, parentid, t), 
                    vector3(frames, from2, to2, parentid, t), 
                    seq
                )
            )
        )
        add_axes!(frames, name, id, AXES_CLASSID_INERTIAL, funs, parentid)

    else
        funs = FrameAxesFunctions{O, N}(
            t -> Rotation{O}(
                twovectors_to_dcm(
                    vector3(frames, from1, to1, parentid, t), 
                    vector3(frames, from2, to2, parentid, t), 
                    seq
                )
            ),
            t -> Rotation{O}(
                twovectors_to_δdcm(
                    vector6(frames, from1, to1, parentid, t), 
                    vector6(frames, from2, to2, parentid, t), 
                    seq
                )
            ),
            t -> Rotation{O}(
                twovectors_to_δ²dcm(
                    vector9(frames, from1, to1, parentid, t), 
                    vector9(frames, from2, to2, parentid, t), 
                    seq
                )
            ),
            t -> Rotation{O}(
                twovectors_to_δ³dcm(
                    vector12(frames, from1, to1, parentid, t), 
                    vector12(frames, from2, to2, parentid, t), 
                    seq
                )
            )
        )
        add_axes!(frames, name, id, AXES_CLASSID_ROTATING, funs, parentid)

    end

end

function _two_vectors_basis(a::AbstractVector, b::AbstractVector, 
    seq::Symbol, fc::Function, fk::Function)

    if seq == :XY
    w = fc(a, b)
    v = fc(w, a)
    u = fk(a)

    elseif seq == :YX
    w = fc(b, a)
    u = fc(a, w)
    v = fk(a)

    elseif seq == :XZ
    v = fc(b, a)
    w = fc(a, v)
    u = fk(a)

    elseif seq == :ZX
    v = fc(a, b)
    u = fc(v, a)
    w = fk(a)

    elseif seq == :YZ
    u = fc(a, b)
    w = fc(u, a)
    v = fk(a)

    elseif seq == :ZY
    u = fc(b, a)
    v = fc(a, u)
    w = fk(a)
    else
    throw(ArgumentError("Invalid rotation sequence."))
    end

    return u, v, w
end

function twovectors_to_dcm(a::AbstractVector, b::AbstractVector, seq::Symbol)
    ut, vt, wt = _two_vectors_basis(a, b, seq, cross3, _keep3)
    u, v, w = unitvec(ut), unitvec(vt), unitvec(wt)

    @inbounds dcm = DCM((u[1], v[1], w[1], u[2], v[2], w[2], u[3], v[3], w[3]))
    return dcm
end

function twovectors_to_δdcm(a::AbstractVector, b::AbstractVector, seq::Symbol)
    u, v, w = _two_vectors_basis(a, b, seq, cross6, _keep6)

    @inbounds @fastmath begin
        ru = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
        rv = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
        rw = sqrt(w[1]^2 + w[2]^2 + w[3]^2)
    end

    δu = δunitvec(u)
    δv = δunitvec(v)
    δw = δunitvec(w)

    @inbounds begin
        dcm = DCM((
            u[1] / ru,
            v[1] / rv,
            w[1] / rw,
            u[2] / ru,
            v[2] / rv,
            w[2] / rw,
            u[3] / ru,
            v[3] / rv,
            w[3] / rw,
        ))

        δdcm = DCM((δu[1], δv[1], δw[1], δu[2], δv[2], δw[2], δu[3], δv[3], δw[3]))
    end

    return dcm, δdcm
end

function twovectors_to_δ²dcm(a::AbstractVector, b::AbstractVector, seq::Symbol)
    u, v, w = _two_vectors_basis(a, b, seq, cross9, _keep9)

    @inbounds @fastmath begin
        ru = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
        rv = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
        rw = sqrt(w[1]^2 + w[2]^2 + w[3]^2)
    end

    δu, δ²u = δunitvec(u), δ²unitvec(u)
    δv, δ²v = δunitvec(v), δ²unitvec(v)
    δw, δ²w = δunitvec(w), δ²unitvec(w)

    @inbounds begin
        dcm = DCM((
            u[1] / ru,
            v[1] / rv,
            w[1] / rw,
            u[2] / ru,
            v[2] / rv,
            w[2] / rw,
            u[3] / ru,
            v[3] / rv,
            w[3] / rw,
        ))

        δdcm = DCM((δu[1], δv[1], δw[1], δu[2], δv[2], δw[2], δu[3], δv[3], δw[3]))

        δ²dcm = DCM((
            δ²u[1], δ²v[1], δ²w[1], δ²u[2], δ²v[2], δ²w[2], δ²u[3], δ²v[3], δ²w[3]
        ))
    end

    return dcm, δdcm, δ²dcm
end

function twovectors_to_δ³dcm(a::AbstractVector, b::AbstractVector, seq::Symbol)
    u, v, w = _two_vectors_basis(a, b, seq, cross12, _keep12)

    @inbounds @fastmath begin
        ru = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
        rv = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
        rw = sqrt(w[1]^2 + w[2]^2 + w[3]^2)
    end

    δu, δ²u, δ³u = δunitvec(u), δ²unitvec(u), δ³unitvec(u)
    δv, δ²v, δ³v = δunitvec(v), δ²unitvec(v), δ³unitvec(v)
    δw, δ²w, δ³w = δunitvec(w), δ²unitvec(w), δ³unitvec(w)

    @inbounds begin
        dcm = DCM((
            u[1] / ru,
            v[1] / rv,
            w[1] / rw,
            u[2] / ru,
            v[2] / rv,
            w[2] / rw,
            u[3] / ru,
            v[3] / rv,
            w[3] / rw,
        ))

        δdcm = DCM((δu[1], δv[1], δw[1], δu[2], δv[2], δw[2], δu[3], δv[3], δw[3]))

        δ²dcm = DCM((
            δ²u[1], δ²v[1], δ²w[1], δ²u[2], δ²v[2], δ²w[2], δ²u[3], δ²v[3], δ²w[3]
        ))

        δ³dcm = DCM((
            δ³u[1], δ³v[1], δ³w[1], δ³u[2], δ³v[2], δ³w[2], δ³u[3], δ³v[3], δ³w[3]
        ))
    end

    return dcm, δdcm, δ²dcm, δ³dcm
end

# These functions are used to guarantee type-stability in the above
for (j, name) in enumerate((:_keep3, :_keep6, :_keep9, :_keep12))
    
    _expr = Expr(:ref, :SA, [Expr(:ref, :a, i) for i in 1:3j]...)

    @eval begin 
        @inline @inbounds function ($name)(a) 
            $_expr
        end
    end

end

