using FrameTransformations
using ReferenceFrameRotations
using StaticArrays
using ForwardDiff
using FiniteDiff
using ChainRulesCore
using Zygote
using Mooncake
using Test

# Build a FrameSystem with analytical point and rotating axes for testing
function build_test_framesystem()
    fr = FrameSystem{2,Float64}()
    add_axes!(fr, :ICRF, 1)

    # Simple rotating axes around Z
    ω = 0.001
    add_axes_rotating!(fr, :RotZ, 2, 1, t -> angle_to_dcm(ω * t, :Z))

    # Root point + dynamical point with known analytical position/velocity
    add_point!(fr, :Root, 1, 1)
    add_point_dynamical!(fr, :Target, 2, 1, 1, t -> SVector(cos(t), sin(t), 1.0))

    return fr
end

const _fr = build_test_framesystem()

# --------------------------------------------------------------------------
# vector3
# --------------------------------------------------------------------------
@testset "vector3 — Zygote vs ForwardDiff" begin
    for t0 in (0.0, 1.0, π / 3)
        fd = ForwardDiff.derivative(t -> sum(vector3(_fr, 1, 2, 1, t)), t0)
        zy = Zygote.gradient(t -> sum(vector3(_fr, 1, 2, 1, t)), t0)[1]
        @test fd ≈ zy rtol = 1e-10
    end
end

@testset "vector3 — Mooncake vs ForwardDiff" begin
    for t0 in (0.0, 1.0, π / 3)
        fd = ForwardDiff.derivative(t -> sum(vector3(_fr, 1, 2, 1, t)), t0)
        f = t -> sum(vector3(_fr, 1, 2, 1, t))
        rule = Mooncake.build_rrule(f, t0)
        _, (_, mc) = Mooncake.value_and_gradient!!(rule, f, t0)
        @test fd ≈ mc rtol = 1e-10
    end
end

# --------------------------------------------------------------------------
# vector6
# --------------------------------------------------------------------------
@testset "vector6 — Zygote vs ForwardDiff" begin
    for t0 in (0.0, 1.0, π / 3)
        fd = ForwardDiff.derivative(t -> sum(vector6(_fr, 1, 2, 1, t)), t0)
        zy = Zygote.gradient(t -> sum(vector6(_fr, 1, 2, 1, t)), t0)[1]
        @test fd ≈ zy rtol = 1e-10
    end
end

@testset "vector6 — Mooncake vs ForwardDiff" begin
    for t0 in (0.0, 1.0, π / 3)
        fd = ForwardDiff.derivative(t -> sum(vector6(_fr, 1, 2, 1, t)), t0)
        f = t -> sum(vector6(_fr, 1, 2, 1, t))
        rule = Mooncake.build_rrule(f, t0)
        _, (_, mc) = Mooncake.value_and_gradient!!(rule, f, t0)
        @test fd ≈ mc rtol = 1e-10
    end
end

# --------------------------------------------------------------------------
# rotation3 — access .m[1] (DCM) which is how downstream code uses it
# --------------------------------------------------------------------------
@testset "rotation3 — Zygote vs ForwardDiff" begin
    for t0 in (0.0, 1.0, π / 3)
        fd = ForwardDiff.derivative(t -> sum(rotation3(_fr, 1, 2, t).m[1]), t0)
        zy = Zygote.gradient(t -> sum(rotation3(_fr, 1, 2, t).m[1]), t0)[1]
        @test fd ≈ zy rtol = 1e-10
    end
end

# --------------------------------------------------------------------------
# rrule correctness: pullback matches FiniteDiff
# --------------------------------------------------------------------------
@testset "rrule pullback vs FiniteDiff" begin
    t0 = 1.0

    # vector3
    g_fd = FiniteDiff.finite_difference_derivative(
        t -> sum(vector3(_fr, 1, 2, 1, t)), t0
    )
    g_zy = Zygote.gradient(t -> sum(vector3(_fr, 1, 2, 1, t)), t0)[1]
    @test g_fd ≈ g_zy rtol = 1e-5

    # vector6
    g_fd6 = FiniteDiff.finite_difference_derivative(
        t -> sum(vector6(_fr, 1, 2, 1, t)), t0
    )
    g_zy6 = Zygote.gradient(t -> sum(vector6(_fr, 1, 2, 1, t)), t0)[1]
    @test g_fd6 ≈ g_zy6 rtol = 1e-5

    # rotation3
    g_fd_r = FiniteDiff.finite_difference_derivative(
        t -> sum(rotation3(_fr, 1, 2, t).m[1]), t0
    )
    g_zy_r = Zygote.gradient(t -> sum(rotation3(_fr, 1, 2, t).m[1]), t0)[1]
    @test g_fd_r ≈ g_zy_r rtol = 1e-5
end
