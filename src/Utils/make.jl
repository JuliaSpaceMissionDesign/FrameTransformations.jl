macro make_struct_fromschema(name, parent, schema...)
    fields=[:($(entry.args[1])::$(entry.args[2])) for entry in schema]
    if parent !== :nothing
        esc(quote struct $name <: $parent
            $(fields...)
            end
        end)
    else
        esc(quote struct $name
            $(fields...)
            end
        end)
    end
end