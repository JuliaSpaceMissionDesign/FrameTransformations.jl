"""
    _two_vectors_basis(a, b, seq::Symbol, fc::Function, fk::Function)

Generate a 3D right-handed orthogonal vector basis and/or its time derivatives from the 
vectors `a` and `b`, according to the directions specified in `seq` and the input cross 
function `fc`. `fk` is a function that filters `a` to guarantee type-stability.

The accepted sequence directions are: `:XY`, `:YX`, `:XZ`, `:ZX`, `:YZ`, `:ZY`

The standard basis, its 1st, 2nd-order and 3rd order time derivatives can be computed by 
passing `cross`, `cross6`, `cross9` or `cross12` to `fc`. The returned vectors will have 
a length of 3, 6 or 9, respectively.

"""
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

"""
    _twovectors_to_dcm(a, b, seq::Symbol, fc::Function, fn::Function, fk::Function)

Generate a direction cosine matrix and/or its time derivatives from the vectors `a` and `b`, 
according to the directions specifeid in `seq`. 

### Notes
`fc` and `fn` are used to control the derivative order. 

"""
function _twovectors_to_dcm(
    a::AbstractVector, b::AbstractVector, seq::Symbol, fc::Function, fn::Function, fk::Function
)
    ut, vt, wt = _two_vectors_basis(a, b, seq, fc, fk)
    u, v, w = fn(ut), fn(vt), fn(wt)

    @inbounds DCM((u[1], v[1], w[1], u[2], v[2], w[2], u[3], v[3], w[3]))
end

"""
    twovectors_to_dcm(a, b, seq)

Generate a direction cosine matrix from two time-dependent vectors `a` and `b`, 
following the directions specified in `seq`. 

### Inputs 
- `a` -- The primary vector that will be aligned with the first directions specified in `seq`. 

- `b` -- The secondary vector. The component of this vector that is orthogonal to the 
         primary vector is aligned with the second direction specified in the sequence `seq`.

- `seq` -- Accepted sequence directions are: 
       `:XY`, `:YX`, `:XZ`, `:ZX`, `:YZ`, `:ZY`


!!! warning 
    The primary and secondary vectors do not have to be orthogonal. However, a great loss of
    precision happens when the two vectors are almost aligned. This function does not perform 
    any check on the angular separation of the two vectors. The user should ensure that the 
    primary and secondary vector differ of at least 1 milliradian.
"""
twovectors_to_dcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross3, unitvec, _keep3)

"""
    twovectors_to_δdcm(a, b, seq)

Compute the time derivative of a direction cosine matrix generated from two time-dependent 
state vectors `a` and `b`, following the directions specified in `seq`. 

### Inputs 
- `a` and `b` -- 6-elements state vectors (position and velocity).
- `seq` -- Accepted sequence directions are: 
       `:XY`, `:YX`, `:XZ`, `:ZX`, `:YZ`, `:ZY`
"""
twovectors_to_δdcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross6, δunitvec, _keep6)

"""
    twovectors_to_δ²dcm(a, b, seq)

Compute the 2nd-order time derivative of a direction cosine matrix generated from two 
time-dependent state vectors `a` and `b`, following the directions specified in `seq`. 

### Inputs 
- `a` and `b` -- 9-elements state vectors (position velocity and acceleration).
- `seq` -- Accepted sequence directions are: `:XY`, `:YX`, `:XZ`, `:ZX`, `:YZ`, `:ZY`
"""
twovectors_to_δ²dcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross9, δ²unitvec, _keep9)

"""
    twovectors_to_δ³dcm(a, b, seq)

Compute the 3rd-order time derivative of a direction cosine matrix generated from two 
time-dependent state vectors `a` and `b`, following the directions specified in `seq`. 

### Inputs 
- `a` and `b` -- 12-elements state vectors (position, velocity, acceleration and jerk).
- `seq` -- Accepted sequence directions are: `:XY`, `:YX`, `:XZ`, `:ZX`, `:YZ`, `:ZY`
"""
twovectors_to_δ³dcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross12, δ³unitvec, _keep12)

"""
    _two_vectors_to_rot6(a, b, seq::Symbol)

Generate a direction cosine matrix and time-derivative, minimising 
the number of repeated computations. 

### See also 
See `twovectors_to_dcm` and `twovectors_to_δdcm` for more information. 
"""
function _two_vectors_to_rot6(
    a::AbstractVector{T}, b::AbstractVector{T}, seq::Symbol
) where {T}
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

"""
    _two_vectors_to_rot9(a, b, seq::Symbol)

Generate a direction cosine matrix and its 1st and 2nd-order time-derivatives, minimising 
the number of repeated computations. 

### See also 
See `twovectors_to_dcm` and `twovectors_to_δ²dcm` for more information. 
"""
function _two_vectors_to_rot9(a::AbstractVector, b::AbstractVector, seq::Symbol)
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

"""
    _two_vectors_to_rot12(a, b, seq::Symbol)

Generate a direction cosine matrix and its 1st, 2nd and 3r-order time-derivatives, minimising 
the number of repeated computations. 

### See also 
See [`twovectors_to_dcm`](@ref) and [`twovectors_to_δdcm`](@ref) for more information. 
"""
function _two_vectors_to_rot12(a::AbstractVector, b::AbstractVector, seq::Symbol)
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
    
    expr = Expr(:ref, :SA, [Expr(:ref, :a, i) for i in 1:3j]...)

    @eval begin 
        @inline @inbounds function ($name)(a) 
            $expr
        end
    end

end
