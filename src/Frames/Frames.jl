module Frames 
    using Logging
    import ForwardDiff.derivative
    using Basic.Utils: format_camelcase
    using Basic.MappedGraphs

    include("rotation.jl")

end