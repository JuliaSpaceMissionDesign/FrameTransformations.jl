using ReferenceFrameRotations
using StaticArrays

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


function vectors_to_xy_dcm(x::AbstractVector, y::AbstractVector)
    uᵤ = x/norm(x)
    w = cross(uᵤ, y)
    wᵤ = w/norm(w)
    v = cross(wᵤ, uᵤ)
    vᵤ = v/norm(v)

    _vectors_to_dcm(uᵤ, vᵤ, wᵤ)
end

function vectors_to_xz_dcm(x::AbstractVector, z::AbstractVector)
    uᵤ = x/norm(x)
    v = cross(z, uᵤ)
    vᵤ = v/norm(v)
    w = cross(uᵤ, vᵤ)
    wᵤ = w/norm(w)

    _vectors_to_dcm(uᵤ, vᵤ, wᵤ)
end

function vectors_to_yz_dcm(y::AbstractVector, z::AbstractVector)
    vᵤ = y/norm(y)
    u = cross(y, z)
    uᵤ = u/norm(u)
    w = cross(uᵤ, vᵤ)
    wᵤ = w/norm(w)

    _vectors_to_dcm(uᵤ, vᵤ, wᵤ)
end

vectors_to_yx_dcm(y, x) = vectors_to_xy_dcm(x, y)
vectors_to_zx_dcm(z, x) = vectors_to_xz_dcm(x, z)
vectors_to_zy_dcm(z, y) = vectors_to_yz_dcm(y, z)

@inbounds function δunit_vector(v::AbstractVector{T}) where T

    x  = v[SA[1, 2, 3]]
    dx = v[SA[4, 5, 6]]

    @fastmath r = norm(x)
    δ = -dot(x, dx)/r^3

    SA{T}[dx[1]/r + δ*x[1], dx[2]/r + δ*x[2], dx[3]/r + δ*x[3]]
end


@fastmath @inbounds function δ²unit_vector(v::AbstractVector{T}) where T

    x  = v[SA[1, 2, 3]]
    dx = v[SA[4, 5, 6]]
    ddx = v[SA[7, 8, 9]]

    r = norm(x)
    dr = norm(dx)

    # δr2 = dx[1]^2+dx[2]^2+dx[3]^2
    # r2 = x[1]^2+x[2]^2+x[3]^2
    # r = sqrt(r2)
    # r3 = r2*r
    # r5 = r3*r2

    δ = dot(x, dx)

    # a = ddx[1]/r - dx[1]/r3*δ - (dx[1]*δ+x[1]*(δr2+dot(x, ddx)))/r3 + 3/r5*x[1]*δ^2
    # b = ddx[2]/r - dx[2]/r3*δ - (dx[2]*δ+x[2]*(δr2+dot(x, ddx)))/r3 + 3/r5*x[2]*δ^2
    # c = ddx[3]/r - dx[3]/r3*δ - (dx[3]*δ+x[3]*(δr2+dot(x, ddx)))/r3 + 3/r5*x[3]*δ^2


    a = ddx[1]/r - dx[1]*dr/r^2 - dx

    SA{T}[a, b, c]
end


@inbounds function cross6(x::AbstractVector, y::AbstractVector)
    u = x[2]*y[3]-x[3]*y[2]
    v = x[3]*y[1]-y[3]*x[1]
    w = x[1]*y[2]-x[2]*y[1]

    δu = x[5]*y[3]+x[2]*y[6]-x[6]*y[2]-x[3]*y[5]
    δv = x[6]*y[1]+x[3]*y[4]-x[4]*y[3]-x[1]*y[6]
    δw = x[4]*y[2]+x[1]*y[5]-x[5]*y[1]-x[2]*y[4]

    SA[u, v, w, δu, δv, δw]
end

@inbounds function cross9(x::AbstractVector, y::AbstractVector)

    u = x[2]*y[3]-x[3]*y[2]
    v = x[3]*y[1]-y[3]*x[1]
    w = x[1]*y[2]-x[2]*y[1]

    δu = x[5]*y[3]+x[2]*y[6]-x[6]*y[2]-x[3]*y[5]
    δv = x[6]*y[1]+x[3]*y[4]-x[4]*y[3]-x[1]*y[6]
    δw = x[4]*y[2]+x[1]*y[5]-x[5]*y[1]-x[2]*y[4]

    δ²u = x[8]*y[3] + 2*x[5]*y[6] + x[2]*y[9] - x[9]*y[2] - 2*x[6]*y[5] - x[3]*y[8]
    δ²v = x[9]*y[1] + 2*x[6]*y[4] + x[3]*y[7] - x[7]*y[3] - 2*x[4]*y[6] - x[1]*y[9]
    δ²w = x[7]*y[2] + 2*x[4]*y[5] + x[1]*y[8] - x[8]*y[1] - 2*x[5]*y[4] - x[2]*y[7]

    SA[u, v, w, δu, δv, δw, δ²u, δ²v, δ²w]
end

function vectors_to_xy_ddcm(x, y)

    w = cross6(x, y);
    v = cross6(w, x);

    δu = δunit_vector(x)
    δv = δunit_vector(v)
    δw = δunit_vector(w)

    _vectors_to_dcm(δu, δv, δw)
end

function vectors_to_xz_ddcm(x, z)

    v = cross6(z, x);
    w = cross6(x, v);

    δu = δunit_vector(x)
    δv = δunit_vector(v)
    δw = δunit_vector(w)

    _vectors_to_dcm(δu, δv, δw)
end

function vectors_to_yz_ddcm(y, z)

    u = cross6(y, z);
    w = cross6(u, y);

    δu = δunit_vector(u)
    δv = δunit_vector(y)
    δw = δunit_vector(w)

    _vectors_to_dcm(δu, δv, δw)
end

vectors_to_yx_ddcm(y, x) = vectors_to_xy_ddcm(x, y)
vectors_to_zx_ddcm(z, x) = vectors_to_xz_ddcm(x, z)
vectors_to_zy_ddcm(z, y) = vectors_to_yz_ddcm(y, z)



function test(t::Number)
    s, c = sincos(t)

    x = SA[c, s, 0]
    y = SA[-s, c, 0]

    vectors_to_zy_dcm(x, y)
end

function dtest(t::Number)
    s, c = sincos(t)

    x = SA[c, s, 0, -s, c, 0]
    y = SA[-s, c, 0, -c, -s, 0]

    vectors_to_zy_ddcm(x, y)
end

function test2(t::Number)
    s, c = sincos(t)

    x = SA[t^2*c, s, 1]
    y = SA[-t*s, c, t^2+3]

    vectors_to_xy_dcm(x, y)
end

function dtest2(t::Number)
    s, c = sincos(t)

    x = SA[t^2*c, s, 1, 2t*c-t^2*s, c, 0]
    y = SA[-t*s, c, t^2+3, -s-t*c, -s, 2t]

    vectors_to_xy_ddcm(x, y)
end

t = π/3
R1 = ForwardDiff.derivative(test2, t)
R2 = dtest2(t)

x = SA[rand(3)...]
y = SA[rand(3)...]
δx = SA[rand(3)...]
δy = SA[rand(3)...]
δδx = SA[rand(3)...]
δδy = SA[rand(3)...]

ux = vcat(x, δx, δδx)
uy = vcat(y, δy, δδy)

A = cross9(ux, uy)[7:9]
B = cross(δδx, y) + 2cross(δx, δy) + cross(x, δδy)
@benchmark cross9($ux, $uy)

function testu(t)
    s, c = sincos(t)
    x = SA[t^2*c, s, 1]
    x/norm(x)
end

function testv(t)
    s, c = sincos(t)
    x = SA[t^2*c, s, 1]
end

t = rand()
s, c = sincos(t)

ForwardDiff.derivative(testu, t)
ForwardDiff.derivative(τ->ForwardDiff.derivative(testu, τ), t)

# y1 = [t^2*c, s, 1, 2t*c-t^2*s, c, 0, 2c-2t*s-t^2*c-2t*s, -s, 0]
y = vcat(testv(t), ForwardDiff.derivative(testv, t),
         ForwardDiff.derivative(τ->ForwardDiff.derivative(testv, τ), t))

δ²unit_vector(y)
ForwardDiff.derivative(τ->ForwardDiff.derivative(testu, τ), t)

# @benchmark dtest($t)
# @benchmark vectors_to_xy_ddcm($x, $y)

# s, c = sincos(t)
# x = [c, s, 0]
# dx = [-s, c, 0]

# y = [-s, c, 0]
# dy = [-c, -s, 0]

# cross(dx, y) + cross(x, dy)
# cross6(SA[x..., dx...], SA[y..., dy...])


# struct Rotation{S, N, T} <: StaticMatrix{N, N, T}
#     data::NTuple{S, DCM{T}}

#     function Rotation{S, T}(x::Tuple{<:DCM}) where {S, T}
#         new{S, 3*S, T}(x)
#     end
# end

# @inline function Base.getindex(rot::Rotation{1, 3, T}, i::Int) where T
#     return rot.data[1][i]
# end

# @inline function Base.getindex(rot::Rotation{2, 6, T}, i::Int) where T
#     row, col = mod(i, 6), i ÷ 6 + 1
#     if row == 0
#         row = 6
#         col -= 1
#     end

#     # println("$i - $row - $col")

#     @inbounds if col < 4 # right
#         row < 4 && return rot.data[1][i - 3*(col-i)]
#         return rot.data[2][i-3col]

#     else # left
#         row < 4 && return T(0)
#         return rot.data[2][i-9-3col]
#     end

# end

# @inline function Base.getindex(rot::Rotation{3, 9, T}, i::Int) where T
#     row, col = i ÷ 9 + 1, mod(i, 9)

#     @inbounds if col < 4 # right
#         row < 4 && return rot.data[1][i - 3*(col-i)]
#         return rot.data[2][i-3col]

#     else # left
#         row < 4 && return T(0)
#         return rot.data[2][i-9-3col]
#     end

# end

# function Tuple(rot::Rotation)
#     return rot.data
# end

# function similar_type(::Type{A}, ::Type{T}, s::Size{(3, 3)}) where {A<:Rotation, T}
#     return Rotation{1, 3, T}
# end



# Rotation(m::DCM{N}) where N = Rotation((m,))

# function Rotation(m::DCM{N}, dm::DCM{N}) where N
#     Rotation((m, dm))
# end

# function Rotation(m::DCM{N}, dm::DCM{N}, ddm::DCM{N}) where N
#     Rotation((m, dm, ddm))
# end

# function Rotation(m::DCM{N}, ω::AbstractVector) where N
#     dm = ddcm(m, SVector(ω))
#     return Rotation{2, N}((m, dm))
# end
