import LinearAlgebra: cross 

"""
    normalize(v::AbstractVector)

Normalise the vector `v`.
"""
@fastmath @inline normalize(v::AbstractVector) = v/norm(v)


"""
    δnormalize(v::AbstractVector)

Compute the time derivative of a unit vector `v`.
"""
function δnormalize(v::AbstractVector{T}) where T

    @inbounds begin 
        r2 = v[1]^2 + v[2]^2 + v[3]^2 
        @fastmath r = sqrt(r2)
        r3 = r2*r

        δ = -(v[1]*v[4] + v[2]*v[5] + v[3]*v[6])/r3

        return SA{T}[v[4]/r + δ*v[1], 
            v[5]/r + δ*v[2], 
            v[6]/r + δ*v[3]]
    end
end


"""
    δ²normalize(v::AbstractVector)

Compute the 2nd-order time derivative of a unit vector `v`.
"""
function δ²normalize(v::AbstractVector{T}) where T

    @inbounds begin 
        δ = v[1]*v[4] + v[2]*v[5] + v[3]*v[6]

        r2 = v[1]^2+v[2]^2+v[3]^2
        @fastmath r  = sqrt(r2)
        r3 = r2*r
        r5 = r3*r2

        Δ = v[1]*v[7] + v[2]*v[8] + v[3]*v[9] + v[4]^2 + v[5]^2 + v[6]^2

        a = v[7]/r - 2v[4]*δ/r3 - v[1]*(Δ/r3 - 3δ^2/r5)
        b = v[8]/r - 2v[5]*δ/r3 - v[2]*(Δ/r3 - 3δ^2/r5)
        c = v[9]/r - 2v[6]*δ/r3 - v[3]*(Δ/r3 - 3δ^2/r5)
    end

    return SA{T}[a, b, c]
end


"""
    cross6(x::AbstractVector, y::AbstractVector)

Compute the cross product between `x` and `y` and its time derivative. 

### Notes 
`x` and `y` must be 6-elements state vectors, with the last elements of each vector 
representing the time derivatives of the first three.
"""
function cross6(x::AbstractVector, y::AbstractVector)

    @inbounds begin 
        u = x[2]*y[3] - x[3]*y[2]
        v = x[3]*y[1] - y[3]*x[1]
        w = x[1]*y[2] - x[2]*y[1]

        δu = x[5]*y[3] + x[2]*y[6] - x[6]*y[2] - x[3]*y[5]
        δv = x[6]*y[1] + x[3]*y[4] - x[4]*y[3] - x[1]*y[6]
        δw = x[4]*y[2] + x[1]*y[5] - x[5]*y[1] - x[2]*y[4]
    end

    return SA[u, v, w, δu, δv, δw]
end


"""
    cross9(x::AbstractVector, y::AbstractVector)

Compute the cross product between `x` and `y` and its 1st and 2nd-order time derivatives. 

### Notes 
`x` and `y` must be 9-elements state vectors (position, velocity and acceleration)
"""
function cross9(x::AbstractVector, y::AbstractVector)

    @inbounds begin 
        u = x[2]*y[3] - x[3]*y[2]
        v = x[3]*y[1] - y[3]*x[1]
        w = x[1]*y[2] - x[2]*y[1]

        δu = x[5]*y[3] + x[2]*y[6] - x[6]*y[2] - x[3]*y[5]
        δv = x[6]*y[1] + x[3]*y[4] - x[4]*y[3] - x[1]*y[6]
        δw = x[4]*y[2] + x[1]*y[5] - x[5]*y[1] - x[2]*y[4]

        δ²u = x[8]*y[3] + 2x[5]*y[6] + x[2]*y[9] - x[9]*y[2] - 2x[6]*y[5] - x[3]*y[8]
        δ²v = x[9]*y[1] + 2x[6]*y[4] + x[3]*y[7] - x[7]*y[3] - 2x[4]*y[6] - x[1]*y[9]
        δ²w = x[7]*y[2] + 2x[4]*y[5] + x[1]*y[8] - x[8]*y[1] - 2x[5]*y[4] - x[2]*y[7]
    end 

    return SA[u, v, w, δu, δv, δw, δ²u, δ²v, δ²w]
end


"""
    _two_vectors_basis(a, b, seq::Symbol, fc::Function)

Generate a 3D right-handed orthogonal vector basis and/or its time derivatives from the 
vectors `a` and `b`, according to the directions specified in `seq` and the input cross 
function `fc`.

The accepted sequence directions are: `:XY`, `:YX`, `:XZ`, `:ZX`, `:YZ`, `:ZY`

The standard basis, its 1st and 2nd-order time derivatives can be computed by passing 
`cross`, `cross6` or `cross9` to `fc`. The returned vectors will have a length of 3, 6 or 9, 
respectively.

"""
function _two_vectors_basis(a::AbstractVector, b::AbstractVector, 
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
        throw(ArgumentError("Invalid rotation sequence."))
    end

    return u, v, w
end


"""
    _twovectors_to_dcm(a, b, seq::Symbol, fc::Function, fn::Function)

Generate a direction cosine matrix and/or its time derivatives from the vectors `a` and `b`, 
according to the directions specifeid in `seq`. 

### Notes
`fc` and `fn` can be used to control the derivative order. 

"""
function _twovectors_to_dcm(a::AbstractVector, b::AbstractVector, seq::Symbol, 
                        fc::Function, fn::Function)

    u, v, w = _two_vectors_basis(a, b, seq, fc)
    u, v, w = fn(u), fn(v), fn(w)

    @inbounds DCM((u[1], v[1], w[1], 
                u[2], v[2], w[2],  
                u[3], v[3], w[3]))

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


### Notes 
The primary and secondary vectors do not have to be orthogonal. However, a great loss of
precision happens when the two vectors are almost aligned. This function does not perform 
any check on the angular separation of the two vectors. The user should ensure that the 
primary and secondary vector differ of at least 1 milliradian.
"""
twovectors_to_dcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross, normalize)


"""
    twovectors_to_ddcm(a, b, seq)

Compute the time derivative of a direction cosine matrix generated from two time-dependent 
state vectors `a` and `b`, following the directions specified in `seq`. 

### Inputs 
- `a` and `b` -- 6-elements state vectors (position and velocity).
- `seq` -- Accepted sequence directions are: 
       `:XY`, `:YX`, `:XZ`, `:ZX`, `:YZ`, `:ZY`
"""
twovectors_to_ddcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross6, δnormalize)

"""
    twovectors_to_ddcm(a, b, seq)

Compute the 2nd-order time derivative of a direction cosine matrix generated from two 
time-dependent state vectors `a` and `b`, following the directions specified in `seq`. 

### Inputs 
- `a` and `b` -- 9-elements state vectors (position and velocity).
- `seq` -- Accepted sequence directions are: `:XY`, `:YX`, `:XZ`, `:ZX`, `:YZ`, `:ZY`
"""
twovectors_to_dddcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross9, δ²normalize)

# Generate a dcm and two identity rotations 
function _two_vectors_to_rot3(a::AbstractVector{T}, b::AbstractVector{T}, seq::Symbol) where T 
    return twovectors_to_dcm(a, b, seq), DCM(T(1)I), DCM(T(1)I)
end 

# Generate a dcm, its derivative and a null rotation
function _two_vectors_to_rot6(a::AbstractVector{T}, b::AbstractVector{T}, seq::Symbol) where T

    u, v, w = _two_vectors_basis(a, b, seq, cross6)

    @inbounds @fastmath begin 
        ru = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
        rv = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
        rw = sqrt(w[1]^2 + w[2]^2 + w[3]^2)
    end

    δu = δnormalize(u)
    δv = δnormalize(v)
    δw = δnormalize(w)

    @inbounds begin 
        dcm = DCM((u[1]/ru, v[1]/rv, w[1]/rw, 
                u[2]/ru, v[2]/rv, w[2]/rw,  
                u[3]/ru, v[3]/rv, w[3]/rw))

        δdcm = DCM((δu[1], δv[1], δw[1], 
                    δu[2], δv[2], δw[2], 
                    δu[3], δv[3], δw[3]))

    end

    return dcm, δdcm, DCM(T(1)I)
end


"""
    _two_vectors_to_rot9(a, b, seq::Symbol)

Generate a direction cosine matrix and its 1st and 2nd-order time-derivatives, minimising 
the number of repeated computations. 

### See also 
See `twovectors_to_dcm` and `twovectors_to_ddcm` for more information. 
"""
function _two_vectors_to_rot9(a::AbstractVector, b::AbstractVector, seq::Symbol)

    u, v, w = _two_vectors_basis(a, b, seq, cross9)

    @inbounds @fastmath begin 
        ru = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
        rv = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
        rw = sqrt(w[1]^2 + w[2]^2 + w[3]^2)
    end

    δu = δnormalize(u)
    δv = δnormalize(v)
    δw = δnormalize(w)

    δ²u = δ²normalize(u)
    δ²v = δ²normalize(v)
    δ²w = δ²normalize(w)

    @inbounds begin 
        dcm = DCM((u[1]/ru, v[1]/rv, w[1]/rw, 
                u[2]/ru, v[2]/rv, w[2]/rw,  
                u[3]/ru, v[3]/rv, w[3]/rw))

        δdcm = DCM((δu[1], δv[1], δw[1], 
                    δu[2], δv[2], δw[2], 
                    δu[3], δv[3], δw[3]))

        δ²dcm = DCM((δ²u[1], δ²v[1], δ²w[1], 
                    δ²u[2], δ²v[2], δ²w[2], 
                    δ²u[3], δ²v[3], δ²w[3]))
    end

    return dcm, δdcm, δ²dcm
end
