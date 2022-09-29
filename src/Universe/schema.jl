export UniverseSchema,
       isvaliduniverse

using JSONSchema: isvalid, Schema, validate

"""
    UniverseSchema

A JSON Schema for input validation purposes.
"""
const UniverseSchema = Schema(
    join(readlines(
        joinpath(@__DIR__, "..", "..", "res", "schemas", "universe.schema.json"))
        )
    )

"""
    isvaliduniverse(x::AbstractDict)::Bool

Validate the universe configuration file.
"""
isvaliduniverse(x::AbstractDict) = validate(UniverseSchema, x) === nothing