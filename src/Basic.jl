module Basic

    using Reexport
    using Logging

    struct AstronautException <: Exception
        msg::String
    end

    Base.showerror(io::IO, err::AstronautException) = print(io, err.msg)

    include(joinpath("Utils", "Utils.jl"))
    @reexport using .Utils

    include(joinpath("Tempo", "Tempo.jl"))
    @reexport using .Tempo

    include(joinpath("Orient", "Orient.jl"))
    @reexport using .Orient

    include(joinpath("Ephem", "Ephem.jl"))
    @reexport using .Ephem
end
