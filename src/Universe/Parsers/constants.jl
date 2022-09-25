export parse_constants!

using Basic.Utils: TPC, load

function parse_constants!(datadict::D1, 
    configdict::D2) where {D1 <: AbstractDict, D2 <: AbstractDict}
    field = configdict["constants"]
    data = load([
        TPC(file) for file in field["ephemeris"]])
    if haskey(field, "database")
        throw(
            error("[Universe] cannot parse constant/database, method not implemented")
        )
    end
    push!(
        datadict,
        :constants => data
    )
end