import LinearAlgebra: cross 

"""
    normalize(v::AbstractVector)

Normalise the vector `v`.
"""
function normalize(v::AbstractVector{T}) where T
    @inbounds begin 
        @fastmath r = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
        SA{T}[v[1]/r, v[2]/r, v[3]/r]
    end
end

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

        SA{T}[v[4]/r + δ*v[1], 
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
        Δ = v[1]*v[7] + v[2]*v[8] + v[3]*v[9]

        r2 = v[1]^2+v[2]^2+v[3]^2
        @fastmath r  = sqrt(r2)
        r3 = r2*r       
        
        dr² = v[4]^2 + v[5]^2 + v[6]^2

        A = (dr² + Δ - 3δ^2/r2)/r3

        a = v[7]/r - 2v[4]*δ/r3 - v[1]*A
        b = v[8]/r - 2v[5]*δ/r3 - v[2]*A
        c = v[9]/r - 2v[6]*δ/r3 - v[3]*A
    end

    SA{T}[a, b, c]
end


"""
    δ³normalize(v::AbstractVector)

Compute the 3rd-order time derivative of a unit vector `v`.
"""
function δ³normalize(v::AbstractVector{T}) where T

    @inbounds begin 
        δ = v[1]*v[4] + v[2]*v[5] + v[3]*v[6]
        Δ = v[1]*v[7] + v[2]*v[8] + v[3]*v[9]
        θ = v[4]*v[7] + v[5]*v[8] + v[6]*v[9]
        φ = v[1]*v[10] + v[2]*v[11] + v[3]*v[12]

        r2 = v[1]^2 + v[2]^2 + v[3]^2
        @fastmath r = sqrt(r2)
        r3 = r2*r

        dr² = v[4]^2 + v[5]^2 + v[6]^2
        δ² = δ^2
        
        A = 3δ/r3
        B = 3*(Δ + dr² - 3δ²/r2)/r3
        C = (3θ + φ - 3δ/r2*(3Δ + 3dr² - 5δ²/r2))/r3

        a = v[10]/r - v[7]*A - v[4]*B - v[1]*C
        b = v[11]/r - v[8]*A - v[5]*B - v[2]*C
        c = v[12]/r - v[9]*A - v[6]*B - v[3]*C

    end

    SA{T}[a, b, c]
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

    SA[u, v, w, δu, δv, δw]
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

    SA[u, v, w, δu, δv, δw, δ²u, δ²v, δ²w]
end


"""
    cross12(x::AbstractVector, y::AbstractVector)

Compute the cross product between `x` and `y` and its 1st, 2nd and 3rd order time derivatives. 

### Notes 
`x` and `y` must be 12-elements state vectors (position, velocity and acceleration)
"""
function cross12(x::AbstractVector, y::AbstractVector)

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

        δ³u =   x[11]*y[3] + 3x[8]*y[6] + 3x[5]*y[9] + x[2]*y[12] +
              - x[12]*y[2] - 3x[9]*y[5] - 3x[6]*y[8] - x[3]*y[11]

        δ³v =   x[12]*y[1] + 3x[9]*y[4] + 3x[6]*y[7] + x[3]*y[10] +
              - x[10]*y[3] - 3x[7]*y[6] - 3x[4]*y[9] - x[1]*y[12]
        
        δ³w =   x[10]*y[2] + 3x[7]*y[5] + 3x[4]*y[8] + x[1]*y[11] +
              - x[11]*y[1] - 3x[8]*y[4] - 3x[5]*y[7] - x[2]*y[10]
    end 

    SA[u, v, w, δu, δv, δw, δ²u, δ²v, δ²w, δ³u, δ³v, δ³w]
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
    twovectors_to_δdcm(a, b, seq)

Compute the time derivative of a direction cosine matrix generated from two time-dependent 
state vectors `a` and `b`, following the directions specified in `seq`. 

### Inputs 
- `a` and `b` -- 6-elements state vectors (position and velocity).
- `seq` -- Accepted sequence directions are: 
       `:XY`, `:YX`, `:XZ`, `:ZX`, `:YZ`, `:ZY`
"""
twovectors_to_δdcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross6, δnormalize)


"""
    twovectors_to_δ²dcm(a, b, seq)

Compute the 2nd-order time derivative of a direction cosine matrix generated from two 
time-dependent state vectors `a` and `b`, following the directions specified in `seq`. 

### Inputs 
- `a` and `b` -- 9-elements state vectors (position velocity and acceleration).
- `seq` -- Accepted sequence directions are: `:XY`, `:YX`, `:XZ`, `:ZX`, `:YZ`, `:ZY`
"""
twovectors_to_δ²dcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross9, δ²normalize)


"""
    twovectors_to_δ³dcm(a, b, seq)

Compute the 3rd-order time derivative of a direction cosine matrix generated from two 
time-dependent state vectors `a` and `b`, following the directions specified in `seq`. 

### Inputs 
- `a` and `b` -- 12-elements state vectors (position, velocity, acceleration and jerk).
- `seq` -- Accepted sequence directions are: `:XY`, `:YX`, `:XZ`, `:ZX`, `:YZ`, `:ZY`
"""
twovectors_to_δ³dcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross12, δ³normalize)


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

    return dcm, δdcm
end


"""
    _two_vectors_to_rot9(a, b, seq::Symbol)

Generate a direction cosine matrix and its 1st and 2nd-order time-derivatives, minimising 
the number of repeated computations. 

### See also 
See `twovectors_to_dcm` and `twovectors_to_δdcm` for more information. 
"""
function _two_vectors_to_rot9(a::AbstractVector, b::AbstractVector, seq::Symbol)

    u, v, w = _two_vectors_basis(a, b, seq, cross9)

    @inbounds @fastmath begin 
        ru = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
        rv = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
        rw = sqrt(w[1]^2 + w[2]^2 + w[3]^2)
    end

    δu, δ²u = δnormalize(u), δ²normalize(u)
    δv, δ²v = δnormalize(v), δ²normalize(v)
    δw, δ²w = δnormalize(w), δ²normalize(w)

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


"""
    _two_vectors_to_rot12(a, b, seq::Symbol)

Generate a direction cosine matrix and its 1st, 2nd and 3r-order time-derivatives, minimising 
the number of repeated computations. 

### See also 
See `twovectors_to_dcm` and `twovectors_to_δdcm` for more information. 
"""
function _two_vectors_to_rot12(a::AbstractVector, b::AbstractVector, seq::Symbol)

    u, v, w = _two_vectors_basis(a, b, seq, cross12)

    @inbounds @fastmath begin 
        ru = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
        rv = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
        rw = sqrt(w[1]^2 + w[2]^2 + w[3]^2)
    end

    δu, δ²u, δ³u = δnormalize(u), δ²normalize(u), δ³normalize(u)
    δv, δ²v, δ³v = δnormalize(v), δ²normalize(v), δ³normalize(v)
    δw, δ²w, δ³w = δnormalize(w), δ²normalize(w), δ³normalize(w)

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

        δ³dcm = DCM((δ³u[1], δ³v[1], δ³w[1], 
                     δ³u[2], δ³v[2], δ³w[2], 
                     δ³u[3], δ³v[3], δ³w[3]))
    end

    return dcm, δdcm, δ²dcm, δ³dcm
end
