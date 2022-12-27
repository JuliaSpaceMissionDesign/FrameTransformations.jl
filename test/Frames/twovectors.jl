import ForwardDiff: derivative
import LinearAlgebra: cross, dot


@testset "TwoVectors" verbose=true begin 
    atol = 1e-10

    # Define two non parallel vectors and their 1st, 2nd and 3rd order time derivatives! 
    get_vector(t) = SA[cos(3t), t*sin(t), t^2*cos(t), 
                    - 3sin(3t), sin(t) + t*cos(t), 2t*cos(t) - t^2*sin(t),
                    - 9cos(3t), 2cos(t) - t*sin(t), 2cos(t) - 4t*sin(t) - t^2*cos(t),
                    27sin(3t), - 3sin(t) - t*cos(t), - 6sin(t) - 6t*cos(t) + t^2*sin(t)]

    get_vector2(t) = SA[t^3, t^2, t, 3t^2, 2t, 1, 6t, 2, 0, 6, 0, 0]

    # Function to compute derivatives of normalization 
    δnorm(t) = get_vector(t)[1:3]/norm(get_vector(t)[1:3])

    # Function to compute derivatives of cross product 
    δcross(t) = cross(get_vector(t)[1:3], get_vector2(t)[1:3])

    # Function to compute derivatives of DCMs 
    δdcm(t, seq) = Frames.twovectors_to_dcm(get_vector(t)[1:3], get_vector2(t)[1:3], seq)

    # Normalisation\cross product routines 
    @testset "Derivatives" verbose=true begin 

        # Check vector normalisation
        a, b = [-2., 3., 1., 4], SA[1, 2, 3]
        @test Frames.normalize(b) ≈ b/norm(b)                   atol=atol
        @test Frames.normalize(a) ≈ SA[a[1:3]...]/norm(a[1:3])  atol=atol
        
        θ = rand(); 

        # Unit vector derivatives 
        # -----------------------

        # Compute 1st, 2nd and 3rd order derivatives of unit vector 
        δu  = derivative(δnorm, θ)
        δ²u = derivative(τ->derivative(δnorm, τ), θ)
        δ³u = derivative(τ->derivative(κ->derivative(δnorm, κ), τ), θ)

        v = get_vector(θ)
        @test Frames.δnormalize(v) ≈ δu     atol=atol
        @test Frames.δ²normalize(v) ≈ δ²u   atol=atol
        @test Frames.δ³normalize(v) ≈ δ³u   atol=atol

        # Cross product derivatives 
        # -------------------------  

        a, b = get_vector(θ), get_vector2(θ)

        for (i, fcn) in enumerate([Frames.cross6, Frames.cross9, Frames.cross12])
            cp = fcn(a, b) 
            
            @test cp[1:3] ≈ cross(a[1:3], b[1:3]) atol=atol # Cross product 
            @test cp[4:6] ≈ derivative(δcross, θ) atol=atol # 1st order derivative   

            if i > 1 # 2nd order derivative 
                @test cp[7:9] ≈ derivative(τ->derivative(δcross, τ), θ) atol=atol      
            end

            if i > 2 # 3rd order derivative 
                @test cp[10:12] ≈ derivative(τ->derivative(κ->derivative(δcross, κ), τ), θ) atol=atol      
            end
        end

    end 

    @testset "Sequence Assembly" begin 

        θ = rand() 
        a, b = get_vector(θ)[1:3], get_vector2(θ)[1:3]
        c = cross(a, b) 

        aᵤ, bᵤ, cᵤ = a/norm(a), b/norm(b), c/norm(c)

        # XY 
        R = Frames.twovectors_to_dcm(a, b, :XY) 

        @test R'*[1, 0, 0] ≈ aᵤ             atol=atol
        @test R'*[0, 0, 1] ≈ cᵤ             atol=atol
        @test dot(R'*[0, 1, 0], aᵤ) ≈ 0     atol=atol

        # YX 
        R = Frames.twovectors_to_dcm(a, b, :YX) 

        @test R'*[0, 1, 0] ≈ aᵤ             atol=atol
        @test R'*[0, 0, 1] ≈ -cᵤ            atol=atol
        @test dot(R'*[1, 0, 0], aᵤ) ≈ 0     atol=atol

        # XZ 
        R = Frames.twovectors_to_dcm(a, b, :XZ) 

        @test R'*[1, 0, 0] ≈ aᵤ             atol=atol
        @test R'*[0, 1, 0] ≈ -cᵤ            atol=atol
        @test dot(R'*[0, 1, 0], aᵤ) ≈ 0     atol=atol

        # ZX 
        R = Frames.twovectors_to_dcm(a, b, :ZX) 

        @test R'*[0, 0, 1] ≈ aᵤ             atol=atol
        @test R'*[0, 1, 0] ≈ cᵤ             atol=atol
        @test dot(R'*[1, 0, 0], aᵤ) ≈ 0     atol=atol

        # YZ 
        R = Frames.twovectors_to_dcm(a, b, :YZ) 

        @test R'*[0, 1, 0] ≈ aᵤ             atol=atol
        @test R'*[1, 0, 0] ≈ cᵤ             atol=atol
        @test dot(R'*[0, 0, 1], aᵤ) ≈ 0     atol=atol

        # ZY
        R = Frames.twovectors_to_dcm(a, b, :ZY) 

        @test R'*[0, 0, 1] ≈ aᵤ             atol=atol
        @test R'*[1, 0, 0] ≈ -cᵤ            atol=atol
        @test dot(R'*[0, 1, 0], aᵤ) ≈ 0     atol=atol

        @test typeof(R) == DCM{Float64}
    end

    @testset "DCM Derivatives" begin 

        θ = rand() 
        a, b = get_vector(θ), get_vector2(θ)

        seq = rand([:XY, :YX, :XZ, :ZX, :YZ, :ZY]) 
        dfun(t) = δdcm(t, seq)

        # 1st order DCM derivative
        B = Frames.twovectors_to_δdcm(a, b, seq)
        @test B ≈ derivative(dfun, θ)                                       atol=atol

        # 2nd order DCM derivative
        C = Frames.twovectors_to_δ²dcm(a, b, seq)
        @test C ≈ derivative(τ->derivative(dfun, τ), θ)                     atol=atol

        # 3rd order DCM derivative
        D = Frames.twovectors_to_δ³dcm(a, b, seq)
        @test D ≈ derivative(τ->derivative(κ->derivative(dfun, κ), τ), θ)   atol=atol
    
        A = Frames.twovectors_to_dcm(a, b, seq) 

        # Assembly of Rotation 6
        R, dR = Frames._two_vectors_to_rot6(a, b, seq)
        @test R     ≈ A    atol=atol 
        @test dR    ≈ B    atol=atol

        # Assembly of Rotation 9
        R, dR, ddR = Frames._two_vectors_to_rot9(a, b, seq)
        @test R     ≈ A    atol=atol 
        @test dR    ≈ B    atol=atol
        @test ddR   ≈ C    atol=atol

        # Assembly of Rotation 12
        R, dR, ddR, dddR = Frames._two_vectors_to_rot12(a, b, seq)
        @test R     ≈ A    atol=atol 
        @test dR    ≈ B    atol=atol
        @test ddR   ≈ C    atol=atol
        @test dddR  ≈ D    atol=atol

    end

end;

