module Universe 

    using OrderedCollections: OrderedDict
    using SHA

    using Basic.Utils
    using CodeGen

    include("Parsers/mappings.jl")
    include("Parsers/bodies.jl")
    include("Parsers/constants.jl")
    include("Parsers/ephemeris.jl")
    include("Parsers/connections.jl")
    include("schema.jl")
    include("parse.jl")

end