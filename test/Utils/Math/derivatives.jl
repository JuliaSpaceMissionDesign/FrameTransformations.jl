
@testset "AD Derivatives" verbose=true begin 

    f = x-> sin(x) + 3x^3
    τ = rand()

    @test D¹(f, τ) ≈ cos(τ) + 9τ^2 atol=1e-14
    @test D²(f, τ) ≈ -sin(τ) + 18τ atol=1e-14
    @test D³(f, τ) ≈ -cos(τ) + 18  atol=1e-14
    
end;