module Basic

    using Reexport

    struct AstronautException <: Exception
        msg::String
    end

    Base.showerror(io::IO, err::AstronautException) = print(io, err.msg)

    include(joinpath("Utils", "Utils.jl"))
    @reexport using .Utils

    include(joinpath("Tempo", "Tempo.jl"))
    @reexport using .Tempo

end
