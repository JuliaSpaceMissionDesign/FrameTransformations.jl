using FrameTransformations

using ReferenceFrameRotations
using LinearAlgebra
using StaticArrays
using Test

@test FrameTransformations.dcm_eltype(DCM{Float64}) == Float64
@test FrameTransformations.promote_dcm_eltype(Tuple{DCM{Int64}, DCM{Float64}}) == Float64

dcm = angle_to_dcm(π/3, :Z)
δdcm = DCM(ddcm(dcm, [0, 0, 1]))
δ²dcm = DCM(ddcm(δdcm, [0, 0, 1]))

r1 = Rotation(dcm, δdcm, δ²dcm)
@test order(r1) == 3
@test size(r1) == (3,)
@test length(r1) == 3
@test r1[1] == dcm 
@test r1[2] == δdcm
@test r1[3] == δ²dcm

@test r1 == Rotation((dcm, δdcm, δ²dcm))

r2 = Rotation{1}(dcm, δdcm, δ²dcm)
@test length(r2) == 1
@test r2[1] == dcm

r3 = Rotation{1}((dcm, δdcm, δ²dcm))
@test length(r3) == 1
@test r3[1] == dcm
@test r2 == r3

r4 = Rotation{2}(1.0I)
@test length(r4) == 2
@test r4[1] == I
@test r4[2] == 0I

r5 = Rotation{3, Float64}(I)
@test length(r5) == 3
@test r5[1] == I
@test r5[2] == 0I
@test r5[3] == 0I

r6 = Rotation{1}(r1)
@test length(r6) == 1
@test r6[1] == dcm

r7 = Rotation{2}(r6)
@test length(r7) == 2
@test r7[1] == dcm 
@test r7[2] == 0I

r8 = Rotation{1, Float32}(r5)
@test length(r8) == 1
@test r8[1] == I

r9 = Rotation{1, Float64}(I)
@test Tuple(r9) == (1, 0, 0, 0, 1, 0, 0, 0, 1)
@test SVector(r9) == [1, 0, 0, 0, 1, 0, 0, 0, 1]
@test inv(r9) == r9
@test inv(r6)[1] == adjoint(dcm)

@test_throws DimensionMismatch r6*r7
r10 = r9*r6
@test r10 == r6

@testset "Constructors" begin

    # Default constructor
    R = Rotation((DCM(1I), DCM(0.0I)))
    @test R[1] == DCM(1I)
    @test R.m[2] == DCM(0I)

    @test typeof(R) == Rotation{2,Float64}

    # dcms constructor 
    A = angle_to_dcm(π / 3, :Z)
    B = DCM(0.0I)
    C = DCM(1.0I)
    R = Rotation(A, B, C)
    @test length(R) == 3

    @test typeof(R) == Rotation{3,Float64}
    for (i, M) in enumerate([A, B, C])
        @test R[i] == M
    end

    # filter constructor 
    R = Rotation{2}(A, B, C)
    @test typeof(R) == Rotation{2,Float64}
    @test R[1] == A
    @test R[2] == B
    @test length(R) == 2

    R = Rotation{5}(A, B, C)
    @test typeof(R) == Rotation{5,Float64}
    @test length(R) == 5
    for (i, M) in enumerate([A, B, C, DCM(0.0I), DCM(0.0I)])
        @test R[i] == M
    end

    # UniformScaling
    R = Rotation{1}(1.0I)
    @test R[1] == C
    @test typeof(R) == Rotation{1,Float64}

    R = Rotation{3,Int64}(I)
    @test typeof(R) == Rotation{3,Int64}
    R[1] == DCM(1I)
    R[2] == DCM(0I)
    R[3] == DCM(0I)

    # Rotation conversion
    R = Rotation(A, B, C)
    R2 = Rotation{1}(R)
    @test R2[1] == A

    # DCM + omega 
    R = Rotation(A, SA[0.0, 0.0, 1.0])
    @test DCM(ddcm(A, SA[0.0, 0.0, 1.0])) == R[2]
end

@testset "Operations" begin
    θ = rand(0:(π / 100):(2π))
    A = angle_to_dcm(θ, :Z)

    R = Rotation(A)
    Ri = Rotation(A')
    R₂ = Rotation(A, DCM(0.0I))

    # Compose rotation
    @test_throws DimensionMismatch R * R₂
    @test (R * Ri)[1] ≈ Rotation{1}(1.0I)[1]
    @test FrameTransformations._compose_rotation(R, Ri)[1] ≈ (R * Ri)[1]

    # Apply rotation
    @testset "apply_rotation" begin 
        t = Translation(1., 0., 0., 0., 1., 0., -1., 0., 0.)
        r1 = Rotation{1, Float64}(I)
        @test_throws DimensionMismatch r1 * t
    
        r2 = Rotation{3}(r1)
        @test r2 * t == t
    
        v = SA[1., 0., 0., 0., 1., 0., -1., 0., 0.]
        @test r2 * v == v
        @test_throws DimensionMismatch r2 * (@SVector zeros(12))
    
        dcm = angle_to_dcm(π/3, :Z)
        r = Rotation(dcm, dcm, dcm)
    
        rt = r * t
        v1 = r[1] * t[1]
        v2 = r[2] * t[1] + r[1] * t[2]  
        v3 = r[3] * t[1] + 2 * r[2] * t[2] + r[1] * t[3]
    
        @test rt[1] == v1
        @test rt[2] == v2
        @test rt[3] == v3
    end

    @test FrameTransformations._apply_rotation(R, SA[1.0, 0.0, 0.0]) ≈ R * SA[1.0, 0.0, 0.0]

    # 2nd order rotations
    A = angle_to_dcm(rand(), :Z)
    B = angle_to_dcm(rand(), rand(), rand(), :XYZ)
    C = angle_to_dcm(rand(), rand(), :YX)
    D = angle_to_dcm(rand(), :X)

    RA = Rotation(A, B)
    RC = Rotation(C, D)

    RO = RA * RC

    @test RO[1] ≈ A * C atol = 1e-12
    @test RO[2] ≈ A * D + B * C atol = 1e-12
    @test typeof(RO) == Rotation{2,Float64}
end
