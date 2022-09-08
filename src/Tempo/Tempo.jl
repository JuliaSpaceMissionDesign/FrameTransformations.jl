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
    include("datetime.jl")
    include("epoch.jl")
    include("offset.jl")
    
    export j2000, j2000seconds, j2000centuries

end