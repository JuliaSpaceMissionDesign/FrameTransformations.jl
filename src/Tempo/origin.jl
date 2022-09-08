abstract type AbstractEpochOrigin end

const EPOCH_ORIGIN = (
    :JulianDate,
    :ModifiedJulianDate,
    :JulianDate2000,
    :ModifiedJulianDate2000
)

const EPOCH_ORIGIN_ACRONYMS = (
    :JD,
    :MJD,
    :J2000,
    :MDJ2000
)

const EPOCH_STARTS = (
    "-4712-01-01T12:00",
    "1858-11-17T00:00" ,
    "2000-01-01T12:00",
    "2000-01-01T00:00"
)

const EPOCH_ORIGIN_TO_J2000 = (
    2.451545e6 * 86400.0,
    51544.5 * 86400.0, 
    0.0, 
    0.5 * 86400.0
)

for (name, acr, off, start) in zip(EPOCH_ORIGIN, EPOCH_ORIGIN_ACRONYMS, 
    EPOCH_ORIGIN_TO_J2000, EPOCH_STARTS)
    acro_str = String(acr)
    name_str = String(name)
    @eval begin
        """
            $($name_str)

        A type representing the $($name_str) ($($acro_str)) epoch origin. 

        With this origin, Epoch reference is `$($start)`.
        """
        struct $name <: AbstractEpochOrigin end

        """
            $($acro_str)

        The singleton instance of the [`$($name_str)`](@ref) epoch origin. 

        With this origin, Epoch reference is `$($start)`.
        """
        const $acr = $name()

        export $name, $acr

        Base.show(io::IO, ::$name) = print(io, "$($acro_str)")
        tryparse(::Val{Symbol($acro_str)}) = $acr
        @inline offset(::$name) = $(off)
    end
end