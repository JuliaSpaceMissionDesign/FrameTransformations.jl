using FrameTransformations

G = FrameSystem{4,Float64}()

add_axes!(G, :ICRF, 1)
add_point!(G, :Earth, 399, 1)

using OrdinaryDiffEq
using DiffEqBase

function two_body!(du, u, p, t)
    μ = p[1]
    r = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
    r3 = r^3
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    du[4] = -μ * u[1] / r3
    du[5] = -μ * u[2] / r3
    du[5] = -μ * u[2] / r3
    nothing
end

u0 = [7000.0, 0.0, 0.0, 0.0, 8.0, 0.1]
p = [398600.0,]
tspan = (0, 10 * 86400.0)

prob = ODEProblem(two_body!, u0, tspan, p)

integrator = DiffEqBase.__init(prob, Vern9(), abstol=1e-14, reltol=1e-12)
add_point_dynamical!(G, :SC, -1, 399, 1, t -> @views(integrator.sol(t)[1:3]), t -> integrator.sol(t))

solve!(integrator)

vector6(G, :Earth, :SC, :ICRF, 123456.0)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
