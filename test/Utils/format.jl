@testset "Function format_camelcase" begin
    @test format_camelcase(String, "snake_case_LINE") == "SnakeCaseLine"
    @test format_camelcase(String, "new-line") == "NewLine"
    @test format_camelcase(String, "CamelLine") == "CamelLine"
end

@testset "Function format_snakecase" begin
    @test format_snakecase(String, "snake_line") == "snake_line"
    @test format_snakecase(String, "new-line") == "new_line"
    @test format_snakecase(String, "CamelCaseLine") == "camel_case_line"
end