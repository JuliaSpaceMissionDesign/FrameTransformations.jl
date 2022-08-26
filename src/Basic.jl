module Basic

    using Reexport
    using Logging

    struct AstronautException <: Exception
        msg::String
    end

    Base.showerror(io::IO, err::AstronautException) = print(io, err.msg)

    include(joinpath("Utils", "Utils.jl"))
    @reexport using .Utils

end
