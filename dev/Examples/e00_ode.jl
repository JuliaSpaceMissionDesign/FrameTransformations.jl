# # [Use Case: Custom propagated orbit](@id example_00_ode)
# _This example was generated on DATEOFTODAY._

# Once the general structure of the [`FrameSystem`](@ref) is understood, we can pass to 
# a use case in which we want to build and exploit our frame system to perform computations 
# inserting a custom orbit. This becomes especially crucial in complex cases like Trajectory 
# Optimization and Navigation Analysis.

# In such scenarios, the trajectory is under design, and the trajectory information might 
# not be completely available. Moreover, in these cases derivatives could be required for 
# various quantities such as time, states, and parameters.

# In this context, we shall remember that `FrameTransformations` is able to perform operations,
# including AD, on the frames **independent variable**, e.g. only time.

# A proper orbit representation is essential to avoid perturbation confusion and ensure 
# proper custom orbit handling.

# ### Frame system setup

# First of all, a new [`FrameSystem`](@ref) shall be created. 

using FrameTransformations

G = FrameSystem{4,Float64}()

add_axes!(G, :ICRF, 1)
add_point!(G, :Earth, 399, 1)

#md # !!! note 
#md #     As the orbit considered in this example is custom, there is no requirement 
#md #     to load ephemeris. 

# ### Custom orbit model

# A custom orbit model is then required. In this case a simple two-body problem is considered 
# but more complex models can be used. We interface with [`SciML`](https://github.com/SciML) 
# ODE solvers to compute the trajectory solution.

using OrdinaryDiffEq
using DiffEqBase

# Now we create the problem with a given initial condition, time-span and gravitational 
# parameter.

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

# It is possible to create a new node using directly the solution of the ODE problem or 
# before the problem is solved, using the integrator:

integrator = DiffEqBase.__init(prob, Vern9(), abstol=1e-14, reltol=1e-12)
add_point_dynamical!(G, :SC, -1, 399, 1, t -> @views(integrator.sol(t)[1:3]), t -> integrator.sol(t))

# Now we can compute the trajectory:

solve!(integrator)

# From now on we can use the frame system with the `:SC` custom node:

vector6(G, :Earth, :SC, :ICRF, 123456.0)