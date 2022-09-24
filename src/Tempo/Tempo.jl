module Tempo

    """
        AbstractDateTimeEpoch
    Supertype for datetime/epoch objects.
    """
    abstract type AbstractDateTimeEpoch end

    """
        TimeScale
    All timescales are subtypes of the abstract type `TimeScale`.
    """
    abstract type TimeScale end

    include("parse.jl")
    include("scales.jl")
    include("origin.jl")
    include("convert.jl")
    include("offset.jl")
    include("datetime.jl")
    include("epoch.jl")

end