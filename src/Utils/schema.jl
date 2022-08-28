export UniverseSchema,
       validate

using JSONSchema: validate, Schema

const UniverseSchema = Schema(
    join(readlines(
        joinpath(@__DIR__, "..", "..", "res", "schemas", "universe.schema.json"))
        )
    )