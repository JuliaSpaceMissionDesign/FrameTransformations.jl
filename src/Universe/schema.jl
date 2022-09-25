export UniverseSchema,
       isvaliduniverse

using JSONSchema: isvalid, Schema, validate

const UniverseSchema = Schema(
    join(readlines(
        joinpath(@__DIR__, "..", "..", "res", "schemas", "universe.schema.json"))
        )
    )

isvaliduniverse(x) = validate(UniverseSchema, x) === nothing