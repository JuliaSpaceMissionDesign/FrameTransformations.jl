using BenchmarkTools 

# CR3BP scheme 
include("frames/Frames.jl")

FRAMES = FrameSystem{Float64}();

in2syn(t::T, y) where T = angle_to_dcm(-t, :Z)

function din2syn(t::T, y) where T 
    s, c = sincos(t)
    DCM{T}(-s, -c, 0, c, -s, 0, 0, 0, 0)
end

icrf2meme = angle_to_dcm(π/3, :Z)
@inertial_axes FRAMES IRF 1 type=InertialReferenceFrame
@inertial_axes FRAMES MEME2000 3 dcm = icrf2meme parent=IRF type=MeanEarthMeanEquinoxJ2000 
@rotating_axes FRAMES SRF 2 IRF in2syn din2syn type=SynodicReferenceFrame

μ = 0.0153 

@root_point FRAMES emb 3 SRF type=EarthMoonBarycenterPoint
@fixed_point FRAMES moon 301 emb SRF SA[1-μ, 0., 0.] type=MoonPoint
@fixed_point FRAMES earth 399 emb SRF SA[-μ, 0., 0.] type=EarthPoint

@benchmark Rotation($FRAMES, $IRF, $SRF, $(0.1))

@benchmark Rotation($FRAMES, $IRF, $MEME2000, $(0.1))

get_vector3(FRAMES, moon, earth, IRF, 0.1)
@benchmark get_vector3($FRAMES, $earth, $moon, $IRF, $(0.1))

function test(axes::AstroAxes, t::Number, y)
    axes.fun(t, y)
end

 
nruns = 100
y = Vector{SVector{3, Float64}}(undef, nruns)
Threads.@threads for i = 1:nruns 
    y[i] = get_vector3(FRAMES, earth, moon, IRF, 0.1)
end

in2syn(t::T) where T = in2syn(t, nothing)


using ForwardDiff
ForwardDiff.derivative(in2syn, π/3)

ForwardDiff.derivative(x->in2syn(x, 0.), π/3)

f = x->in2syn(x, 0.)
@benchmark ForwardDiff.derivative($f, π/3)

function in2syn(t::T, y) where T 
    angle_to_dcm(-t, :Z)
end    

function test(t::T, y) where T 
    ForwardDiff.derivative(x->in2syn(x, y), t)
end

test(π/3, stv)

stv = SA[zeros(6)...]
@benchmark test(π/3, $stv)

function test(x)
    x+2
end

ForwardDiff.derivative(test, 2.)

struct Test{T}
    fun::FunctionWrapper{T, Tuple{T}}
end

T = Float64
x = Test{ForwardDiff.Dual{ForwardDiff.Tag{typeof(test), T}, 1}}(test);
x.fun(2.)

x = Test{Any}(test)
ForwardDiff.derivative(x.fun, 2.)
@benchmark ForwardDiff.derivative($x.fun, $2.)


