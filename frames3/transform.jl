using ReferenceFrameRotations
using StaticArrays
using ForwardDiff

struct Rotation{S, N}
    m::NTuple{S, DCM{N}}
end

function Rotation(m::DCM{N}) where N
    Rotation((m,))
end

function Rotation(m::DCM{N}, dm::DCM{N}) where N
    Rotation((m, dm))
end

function Rotation(m::DCM{N}, dm::DCM{N}, ddm::DCM{N}) where N
    Rotation((m, dm, ddm))
end

function Rotation(m::DCM{N}, ω::AbstractVector) where N
    dm = ddcm(m, SVector(ω))
    return Rotation{2, N}((m, dm))
end


function _vectors_to_dcm(x::AbstractVector, y::AbstractVector, z::AbstractVector)
    @inbounds DCM((x[1], y[1], z[1], x[2], y[2], z[2], x[3], y[3], z[3]))
end

unit_vector(v::AbstractVector) = v/norm(v)

@inbounds function δunit_vector(v::AbstractVector{T}) where T

    r2 = v[1]^2 + v[2]^2 + v[3]^2 
    @fastmath r = sqrt(r2)
    r3 = r2*r

    δ = -(v[1]*v[4] + v[2]*v[5] + v[3]*v[6])/r3

    SA{T}[v[4]/r + δ*v[1], 
          v[5]/r + δ*v[2], 
          v[6]/r + δ*v[3]]
end

@inbounds function δ²unit_vector(v::AbstractVector{T}) where T

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

@inbounds function cross6(x::AbstractVector, y::AbstractVector)

    u = x[2]*y[3] - x[3]*y[2]
    v = x[3]*y[1] - y[3]*x[1]
    w = x[1]*y[2] - x[2]*y[1]

    δu = x[5]*y[3] + x[2]*y[6] - x[6]*y[2] - x[3]*y[5]
    δv = x[6]*y[1] + x[3]*y[4] - x[4]*y[3] - x[1]*y[6]
    δw = x[4]*y[2] + x[1]*y[5] - x[5]*y[1] - x[2]*y[4]

    SA[u, v, w, δu, δv, δw]
end


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

function _vec_xy_dcm(x::AbstractVector, y::AbstractVector, 
                           fc::Function, uv::Function)

    w = fc(x, y)
    v = fc(w, x)

    _vectors_to_dcm(uv(x), uv(v), uv(w))

end

function _vec_xz_dcm(x::AbstractVector, z::AbstractVector, 
                     fc::Function, uv::Function)

    v = fc(z, x)
    w = fc(x, v)

    _vectors_to_dcm(uv(x), uv(v), uv(w))

end

function _vec_yz_dcm(y::AbstractVector, z::AbstractVector, 
                     fc::Function, uv::Function)

    u = fc(y, z)
    w = fc(u, y)

    _vectors_to_dcm(uv(u), uv(y), uv(w))

end




function vectors_to_yz_dcm(y::AbstractVector, z::AbstractVector)

    u = cross(y, z)
    w = cross(uᵤ, vᵤ)

    uᵤ = u/norm(u)
    vᵤ = y/norm(y) 
    wᵤ = w/norm(w)

    _vectors_to_dcm(uᵤ, vᵤ, wᵤ)
end

vectors_to_yx_dcm(y, x) = vectors_to_xy_dcm(x, y)
vectors_to_zx_dcm(z, x) = vectors_to_xz_dcm(x, z)
vectors_to_zy_dcm(z, y) = vectors_to_yz_dcm(y, z)

function d_vectors_to_xy_dcm(x, y)

    w = cross6(x, y);
    v = cross6(w, x);

    δu = δunit_vector(x)
    δv = δunit_vector(v)
    δw = δunit_vector(w)

    _vectors_to_dcm(δu, δv, δw)
end

function d_vectors_to_xz_dcm(x, z)

    v = cross6(z, x);
    w = cross6(x, v);

    δu = δunit_vector(x)
    δv = δunit_vector(v)
    δw = δunit_vector(w)

    _vectors_to_dcm(δu, δv, δw)
end

function d_vectors_to_yz_dcm(y, z)

    u = cross6(y, z);
    w = cross6(u, y);

    δu = δunit_vector(u)
    δv = δunit_vector(y)
    δw = δunit_vector(w)

    _vectors_to_dcm(δu, δv, δw)
end

d_vectors_to_yx_dcm(y, x) = d_vectors_to_xy_dcm(x, y)
d_vectors_to_zx_dcm(z, x) = d_vectors_to_xz_dcm(x, z)
d_vectors_to_zy_dcm(z, y) = d_vectors_to_yz_dcm(y, z)


function dd_vectors_to_xy_dcm(x, y)

    w = cross9(x, y);
    v = cross9(w, x);

    δ²u = δ²unit_vector(x)
    δ²v = δ²unit_vector(v)
    δ²w = δ²unit_vector(w)

    _vectors_to_dcm(δ²u, δ²v, δ²w)
end

function dd_vectors_to_xz_ddcm(x, z)

    v = cross9(z, x);
    w = cross9(x, v);

    δ²u = δ²unit_vector(x)
    δ²v = δ²unit_vector(v)
    δ²w = δ²unit_vector(w)

    _vectors_to_dcm(δ²u, δ²v, δ²w)
end

function dd_vectors_to_yz_ddcm(y, z)

    u = cross9(y, z);
    w = cross9(u, y);

    δ²u = δ²unit_vector(u)
    δ²v = δ²unit_vector(y)
    δ²w = δ²unit_vector(w)

    _vectors_to_dcm(δ²u, δ²v, δ²w)
end

dd_vectors_to_yx_dcm(y, x) = dd_vectors_to_xy_dcm(x, y)
dd_vectors_to_zx_dcm(z, x) = dd_vectors_to_xz_dcm(x, z)
dd_vectors_to_zy_dcm(z, y) = dd_vectors_to_yz_dcm(y, z)


function test(t::Number)
    s, c = sincos(t)

    x = SA[t^2*c, s, 1]
    y = SA[-t*s, c, t^2+3]

    vectors_to_xy_dcm(x, y)
end

function d_test(t::Number)
    s, c = sincos(t)

    x = SA[t^2*c, s, 1, 2t*c-t^2*s, c, 0]
    y = SA[-t*s, c, t^2+3, -s-t*c, -s, 2t]

    d_vectors_to_xy_dcm(x, y)
end

function dd_test(t::Number)
    s, c = sincos(t)

    x = SA[t^2*c, s, 1, 2t*c-t^2*s, c, 0, 2c - 4t*s - t^2*c, -s, 0]
    y = SA[-t*s, c, t^2+3, -s-t*c, -s, 2t, -2c+t*s, -c, 2]

    dd_vectors_to_xy_dcm(x, y)
end

# First order derivative test 
t = π/3;
R1 = ForwardDiff.derivative(test, t);
R2 = d_test(t);
maximum(abs.(R1-R2))

# Second order derivative test 
R1 = ForwardDiff.derivative(τ->ForwardDiff.derivative(test, τ), t);
R2 = dd_test(t);
maximum(abs.(R1-R2))

# Zero, First and Second order time derivatives of cross product test 
x, y = SA[rand(3)...], SA[rand(3)...];
δx, δy = SA[rand(3)...], SA[rand(3)...];
δδx, δδy = SA[rand(3)...], SA[rand(3)...];
ux, uy = vcat(x, δx, δδx), vcat(y, δy, δδy);

A = cross9(ux, uy)[7:9];
B = cross(δδx, y) + 2cross(δx, δy) + cross(x, δδy);
maximum(abs.(A-B))

function testvx(t)
    s, c = sincos(t)
    SA[t^2*c, s, 1]
end

function testvy(t)
    s, c = sincos(t)
    SA[-t*s, c, t^2+3]
end

function testux(t)
    s, c = sincos(t)
    x = SA[t^2*c, s, 1]
    x/norm(x)
end

function testuy(t)
    s, c = sincos(t)
    x = SA[-t*s, c, t^2+3]
    x/norm(x)
end

t = rand()
s, c = sincos(t)

x1 = [t^2*c, s, 1, 2t*c-t^2*s, c, 0, 2c - 4t*s - t^2*c, -s, 0];
x2 = vcat(testvx(t), ForwardDiff.derivative(testvx, t),
         ForwardDiff.derivative(τ->ForwardDiff.derivative(testvx, τ), t));

maximum(abs.(x1-x2))

y1 = [-t*s, c, t^2+3, -s-t*c, -s, 2t, -2c+t*s, -c, 2];
y2 = vcat(testvy(t), ForwardDiff.derivative(testvy, t),
         ForwardDiff.derivative(τ->ForwardDiff.derivative(testvy, τ), t));

maximum(abs.(y1-y2))

u1 = δ²unit_vector(x1);
u2 = ForwardDiff.derivative(τ->ForwardDiff.derivative(testux, τ), t);
maximum(abs.(u1-u2))

v1 = δ²unit_vector(y1);
v2 = ForwardDiff.derivative(τ->ForwardDiff.derivative(testuy, τ), t);
maximum(abs.(v1-v2))

# δ²x/r + 3x*(x*δx + y*δy + z*δz)^2/(r^5) - (x*(2δx^2 +2δy^2 + 2δz^2 + 4x*δx^2))

function testv(λ::Number)

    s, c = sincos(λ)
    

end