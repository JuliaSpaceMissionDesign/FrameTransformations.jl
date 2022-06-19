module Tempo

    """
        AbstractDateTimeEpoch
    Supertype for datetime/epoch objects
    """
    abstract type AbstractDateTimeEpoch end

    """
        TimeScale
    All timescales are subtypes of the abstract type `TimeScale`.
    """
    abstract type TimeScale end

    include("instant.jl")
    include("parse.jl")
    include("datetime.jl")
    include("epoch.jl")
    include("scales.jl")

end