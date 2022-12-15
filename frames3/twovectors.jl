

# compute unit vector 
@fastmath @inline normalize(v::AbstractVector) = v/norm(v)

# comput time-derivative of unit vector 
@inbounds function δnormalize(v::AbstractVector{T}) where T

    r2 = v[1]^2 + v[2]^2 + v[3]^2 
    @fastmath r = sqrt(r2)
    r3 = r2*r

    δ = -(v[1]*v[4] + v[2]*v[5] + v[3]*v[6])/r3

    SA{T}[v[4]/r + δ*v[1], 
          v[5]/r + δ*v[2], 
          v[6]/r + δ*v[3]]
end

# compute 2nd order time derivative of unit vector 
@inbounds function δ²normalize(v::AbstractVector{T}) where T

    δ = v[1]*v[4] + v[2]*v[5] + v[3]*v[6]

    r2 = v[1]^2+v[2]^2+v[3]^2
    @fastmath r  = sqrt(r2)
    r3 = r2*r
    r5 = r3*r2

    Δ = v[1]*v[7] + v[2]*v[8] + v[3]*v[9] + v[4]^2 + v[5]^2 + v[6]^2

    a = v[7]/r - 2v[4]*δ/r3 - v[1]*Δ/r3 + 3v[1]*δ^2/r5
    b = v[8]/r - 2v[5]*δ/r3 - v[2]*Δ/r3 + 3v[2]*δ^2/r5
    c = v[9]/r - 2v[6]*δ/r3 - v[3]*Δ/r3 + 3v[3]*δ^2/r5

    SA{T}[a, b, c]
end

# cross-product and its derivative 
@inbounds function cross6(x::AbstractVector, y::AbstractVector)

    u = x[2]*y[3] - x[3]*y[2]
    v = x[3]*y[1] - y[3]*x[1]
    w = x[1]*y[2] - x[2]*y[1]

    δu = x[5]*y[3] + x[2]*y[6] - x[6]*y[2] - x[3]*y[5]
    δv = x[6]*y[1] + x[3]*y[4] - x[4]*y[3] - x[1]*y[6]
    δw = x[4]*y[2] + x[1]*y[5] - x[5]*y[1] - x[2]*y[4]

    SA[u, v, w, δu, δv, δw]
end

# Zero, First and Second order time derivatives of cross product test 
@inbounds function cross9(x::AbstractVector, y::AbstractVector)

    u = x[2]*y[3] - x[3]*y[2]
    v = x[3]*y[1] - y[3]*x[1]
    w = x[1]*y[2] - x[2]*y[1]

    δu = x[5]*y[3] + x[2]*y[6] - x[6]*y[2] - x[3]*y[5]
    δv = x[6]*y[1] + x[3]*y[4] - x[4]*y[3] - x[1]*y[6]
    δw = x[4]*y[2] + x[1]*y[5] - x[5]*y[1] - x[2]*y[4]

    δ²u = x[8]*y[3] + 2x[5]*y[6] + x[2]*y[9] - x[9]*y[2] - 2x[6]*y[5] - x[3]*y[8]
    δ²v = x[9]*y[1] + 2x[6]*y[4] + x[3]*y[7] - x[7]*y[3] - 2x[4]*y[6] - x[1]*y[9]
    δ²w = x[7]*y[2] + 2x[4]*y[5] + x[1]*y[8] - x[8]*y[1] - 2x[5]*y[4] - x[2]*y[7]

    SA[u, v, w, δu, δv, δw, δ²u, δ²v, δ²w]
end

function _twovectors_to_dcm(a::AbstractVector, b::AbstractVector, seq::Symbol, 
                            fc::Function, fn::Function)

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

    end 

    u, v, w = fn(u), fn(v),  fn(w)
    @inbounds DCM((u[1], v[1], w[1], u[2], v[2], w[2],  u[3], v[3], w[3]))

end

twovectors_to_dcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross, normalize)
twovectors_to_ddcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross6, δnormalize)
twovectors_to_dddcm(a, b, seq) = _twovectors_to_dcm(a, b, seq, cross9, δ²normalize)