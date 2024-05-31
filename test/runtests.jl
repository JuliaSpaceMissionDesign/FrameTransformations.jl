using Test 

const GROUP = get(ENV, "GROUP", "All")

@time begin
    if GROUP == "All" || GROUP == "Core"
        include("Core/Core.jl")
    end
    if GROUP == "All" || GROUP == "Definitions"
        include("Definitions/Definitions.jl")
    end
end