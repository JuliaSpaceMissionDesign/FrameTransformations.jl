export AstronautGenericError, NotImplementedError

"""
    AstronautGenericException

A supertype for all astronaut related errors.
"""
abstract type AstronautGenericException <: Exception end

"""
    @create_module_error

Create a type representing an error.
"""
macro create_module_error(name, supertype, descr)
    return esc(
        quote
            """
                $($(name))

            A type representing $($(descr))
            """
            struct $name <: $supertype
                modl::String
                msg::String
            end
            function Base.showerror(io::IO, err::$name)
                return print(io, "[", err.modl, "] ", err.msg)
            end
        end,
    )
end

#
# error types 
# 

@create_module_error AstronautGenericError AstronautGenericException "generic errors"
@create_module_error NotImplementedError AstronautGenericException "not implemented errors"
@create_module_error EphemerisError AstronautGenericException "ephemeris errors"
