@testset "Abstract load" begin
    @test load() === nothing
end

@testset "Function load(JSON)" begin
    @test load(JSON(joinpath(@__DIR__, "assets", "test_json.json")))[:check] == 1
end

@testset "Function load(TXT)" begin
    @test load(TXT(joinpath(@__DIR__, "assets", "test_txt.txt")))[1] == "check: 1"
end

@testset "Function load(YAML)" begin
    @test load(YML(joinpath(@__DIR__, "assets", "test_yaml.yml")))[:check] == 1

    @test load(YAML(joinpath(@__DIR__, "assets", "test_yaml.yml")))[:check] == 1
end

@testset "Function load(TPC)" begin
    file = TPC(joinpath(@__DIR__, "..", "assets", "gms.tpc"))
    data = load(file)
    @test data[399][:gm] == 3.9860043543609598E+05
    @test_throws KeyError data[-100]

    mapped = Dict{Int64,Dict{Symbol,Union{Float64,Vector{Float64}}}}()
    Utils.load_tpc!(mapped, Utils.filepath(file))
    @test data == mapped
end

@testset "Function load(Vector(TPC))" begin
    gms = TPC(joinpath(@__DIR__, "..", "assets", "gms.tpc"))
    orient = TPC(joinpath(@__DIR__, "..", "assets", "orient.tpc"))

    data = load([gms, orient])
    for (k, l, fv, lv) in zip(
        [:pole_ra, :pole_dec, :pm, :nut_prec_ra, :nut_prec_dec, :nut_prec_pm],
        [3, 3, 3, 15, 15, 15],
        [268.056595, 64.495303, 284.95, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.002150, 0.000926, 0.0],
    )
        @test length(data[599][k]) == l # test length
        @test data[599][k][1] == fv  # test first value
        @test data[599][k][end] == lv  # test last value
    end
end
