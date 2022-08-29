export parse_mappings, parse_mappings!

using OrderedCollections: OrderedDict

function parse_mappings(configdict::D) where {D <: AbstractDict{Symbol, Any}}
    parsed = OrderedDict{Symbol, Any}()
    for line in configdict[:mappings] 
        parser = Symbol("parse_$(line[:name])")
        @eval begin
            push!(
                $parsed,
                Symbol($(line[:name])) => $parser($(line[:path]))
            )
        end
    end
    parsed
end

function parse_mappings!(datadict::D, configdict::D) where {D <: AbstractDict{Symbol, Any}}
    parsed = parse_mappings(configdict)
    push!(
        datadict, 
        :mappings => parsed
    )
end