
@testset "Rotation" verbose = true begin
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
        @test Frames._compose_rot(R, Ri)[1] ≈ (R * Ri)[1]

        # Apply rotation
        @test Frames._apply_rot(R, [1.0, 0.0, 0.0]) ≈ SA[cos(θ), -sin(θ), 0.0]
        @test Frames._apply_rot(R, SA[1.0, 0.0, 0.0]) ≈ R * SA[1.0, 0.0, 0.0]
        @test Frames._apply_rot(R, [1.0, 0.0, 0.0]) ≈ R * [1.0, 0.0, 0.0]

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

    @testset "Utility functions" begin
        R = Rotation((DCM(1I), DCM(0.0I)))
        @test Frames.order(R) == 2
        @test size(R) == (6, 6)
        @test StaticArrays.Size(R) == (6, 6)
        @test R[1] == R.m[1]
        @test typeof(R) == Rotation{2,Float64}

        @test StaticArrays.similar_type(R) == Rotation{2,Float64}
        @test StaticArrays.similar_type(R, BigFloat) == Rotation{2, BigFloat}

        # Convert to tuple
        A = angle_to_dcm(rand(), rand(), rand(), :XYZ)
        @test Tuple(Rotation(A)) === Tuple(A)

        # Convert to SMatrix
        @test typeof(SMatrix(R)) == SMatrix{6,6,Float64,36}

        # Convert to MMatrix
        @test typeof(MMatrix(R)) == MMatrix{6,6,Float64,36}

        # Inverse rotation
        A = angle_to_dcm(π / 3, :Z)
        R = Rotation(A)
        @test inv(R)[1] ≈ A'
        @test Frames._inverse_rot(R) == inv(R)
    end
end;
