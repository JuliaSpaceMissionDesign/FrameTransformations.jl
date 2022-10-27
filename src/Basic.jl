module Basic

    using Reexport
    using Logging

    struct AstronautException <: Exception
        msg::String
    end

    Base.showerror(io::IO, err::AstronautException) = print(io, err.msg)
    
    # Common 
    include("graph.jl")

    include(joinpath("Utils", "Utils.jl"))
    @reexport using .Utils

    include(joinpath("Tempo", "Tempo.jl"))
    @reexport using .Tempo

    include(joinpath("Bodies", "Bodies.jl"))
    @reexport using .Bodies

    include(joinpath("Ephemeris", "Ephemeris.jl"))
    @reexport using .Ephemeris

    include(joinpath("Orient", "Orient.jl"))
    @reexport using .Orient

    # include(joinpath("Frames", "Frames.jl"))
    # @reexport using .Frames

    include(joinpath("Universe", "Universe.jl"))
    @reexport using .Universe

end
