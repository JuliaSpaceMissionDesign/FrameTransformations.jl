
@testset "Vectors" verbose = true begin
    atol = 1e-12

    # Define two non parallel vectors and their 1st, 2nd and 3rd order time derivatives! 
    function get_vector(t)
        return SA[
            cos(3t),
            t * sin(t),
            t^2 * cos(t),
            -3sin(3t),
            sin(t) + t * cos(t),
            2t * cos(t) - t^2 * sin(t),
            -9cos(3t),
            2cos(t) - t * sin(t),
            2cos(t) - 4t * sin(t) - t^2 * cos(t),
            27sin(3t),
            -3sin(t) - t * cos(t),
            -6sin(t) - 6t * cos(t) + t^2 * sin(t),
        ]
    end

    get_vector2(t) = SA[t^3, t^2, t, 3t^2, 2t, 1, 6t, 2, 0, 6, 0, 0]

    # Function to compute derivatives of normalization 
    δnorm(t) = get_vector(t)[1:3] / norm(get_vector(t)[1:3])

    # Function to compute derivatives of cross product 
    δcross(t) = cross(get_vector(t)[1:3], get_vector2(t)[1:3])

    # Normalisation\cross product routines 

    # Unit vector derivatives 
    # -----------------------
    @testset "Normalisation" verbose = true begin

        # Check vector normalisation
        a, b = [-2.0, 3.0, 1.0, 4], SA[1, 2, 3]
        @test Utils.normalize(b) ≈ b / norm(b) atol = atol
        @test Utils.normalize(a) ≈ SA[a[1:3]...] / norm(a[1:3]) atol = atol

        θ = rand()

        # Compute 1st, 2nd and 3rd order derivatives of unit vector 
        v = get_vector(θ)
        @test Utils.δnormalize(v) ≈ D¹(δnorm, θ) atol = atol
        @test Utils.δ²normalize(v) ≈ D²(δnorm, θ) atol = atol
        @test Utils.δ³normalize(v) ≈ D³(δnorm, θ) atol = atol
    end

    # Cross product derivatives 
    # -------------------------  
    @testset "Cross Product" verbose = true begin
        θ = rand()

        a, b = get_vector(θ), get_vector2(θ)

        for (i, fcn) in enumerate([Utils.cross6, Utils.cross9, Utils.cross12])
            cp = fcn(a, b)

            @test cp[1:3] ≈ cross(a[1:3], b[1:3]) atol = atol # Cross product 
            @test cp[4:6] ≈ D¹(δcross, θ) atol = atol # 1st order derivative   

            if i > 1 # 2nd order derivative 
                @test cp[7:9] ≈ D²(δcross, θ) atol = atol
            end

            if i > 2 # 3rd order derivative 
                @test cp[10:12] ≈ D³(δcross, θ) atol = atol
            end
        end
    end
end
