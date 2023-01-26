
@testset "AD Derivatives" verbose=true begin 

    f = x-> sin(x) + 3x^3
    τ = rand()

    @test D¹(f, τ) ≈ cos(τ) + 9τ^2
    @test D²(f, τ) ≈ -sin(τ) + 18τ
    @test D³(f, τ) ≈ -cos(τ) + 18
    
end;