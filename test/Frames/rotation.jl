using StaticArrays
using ReferenceFrameRotations

@testset "Rotation" verbose=true begin
    @testset "Constructors" begin
    
        # Default constructor
        R = Rotation((DCM(1I), DCM(0.0I)))
        @test R[1] == DCM(1I)
        @test R.m[2] == DCM(0I)

        # dcms constructor 
        A = angle_to_dcm(π/3, :Z)
        B = DCM(0.0I)
        C = DCM(1.0I)
        R = Rotation(A, B, C)
        @test R[3] == C

        # UniformScaling
        R = Rotation{1}(1.0I) 
        @test R[1] == C

        R = Rotation{1, Int64}(I)
        @test typeof(R) == Rotation{1, Int64}

        # Rotation 
        R = Rotation(A, B, C)
        R2 = Rotation{1}(R)
        @test R2[1] == A
        @test_throws DimensionMismatch Rotation{5}(R)

        # DCM + omega 
        R = Rotation(A, SA[0.0, 0.0, 1.0])
        @test DCM(ddcm(A, SA[0.0, 0.0, 1.0])) == R[2]
    end

    @testset "Operations" begin
        
        θ = rand(0:π/100:2*π)
        A = angle_to_dcm(θ, :Z)
        R = Rotation(A)
        R1 = Rotation(A')
        R2 = Rotation(A, DCM(0.0I))

        # Compose rotation
        @test_throws DimensionMismatch R * R2
        @test (R * R1)[1] ≈ Rotation{1}(1.0I)[1]
        @test Frames._compose_rot(R, R1)[1] ≈ (R * R1)[1]

        # Apply rotation
        @test Frames._apply_rot(R, [1.0, 0.0, 0.0]) ≈ SA[cos(θ), -sin(θ), 0.0]
        @test Frames._apply_rot(R, SA[1.0, 0.0, 0.0]) ≈ R * SA[1.0, 0.0, 0.0]
        @test Frames._apply_rot(R, [1.0, 0.0, 0.0]) ≈ R * [1.0, 0.0, 0.0]
        
    end
    
    @testset "Utility functions" begin
        
        R = Rotation((DCM(1I), DCM(0.0I)))
        @test Frames.order(R) == 2
        @test size(R) == (6, 6)
        @test getindex(R, 1) == R.m[1]
        @test typeof(R) == Rotation{2, Float64}

        # Convert to tuple 
        R = Rotation(DCM(1.0I))
        @test Tuple(R) === (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

        # Convert to SMatrix
        @test typeof(SMatrix(R)) == SMatrix{3, 3, Float64, 9}

        # Convert to MMatrix
        @test typeof(MMatrix(R)) == MMatrix{3, 3, Float64, 9}
        
        # Inverse rotation
        A = angle_to_dcm(π/3, :Z)
        R = Rotation(A)
        @test inv(R)[1] ≈ angle_to_dcm(π/3, :Z)'
        @test Frames._inverse_rot(R) == inv(R)

    end    
end