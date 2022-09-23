export DJ2000, DJM0, DMJD

abstract type AbstractEpochOrigin end

"""
    DJ2000

Reference epoch (J2000.0), Julian Date (`2451545.0`). 
It is `12:00 01-01-2000`.
"""
const DJ2000 = 2451545.0

"""
    DMJD

Reference epoch (J2000.0), Modified Julian Date (`51544.5`).
"""
const DMJD = 51544.5

"""
    DJM0

Julian Date of Modified Julian Date zero (`2400000.5`).
It is `00:00 17-11-1858`.
"""
const DJM0 = 2400000.5

const EPOCH_ORIGIN = (
    :JulianDate,
    :ModifiedJulianDate,
    :JulianDate2000,
    :ModifiedJulianDate2000,
)

const EPOCH_ORIGIN_ACRONYMS = (
    :JD,
    :MJD,
    :J2000,
    :MJ2000
)

const EPOCH_STARTS = (
    "-4712-01-01T12:00",
    "1858-11-17T00:00" ,
    "2000-01-01T12:00",
    "2000-01-01T00:00"
)

const EPOCH_ORIGIN_TO_J2000 = (
    DJ2000,
    DMJD, 
    0.0, 
    -0.5,
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
        Base.tryparse(::Val{Symbol($acro_str)}) = $acr

        """
            offset(::$($name))

        Offset to shift J2000 epochs (with origin at `2000-01-01T12:00`) 
        to [`$($name_str)`](@ref) (with origin at `$($start)`)
        """
        @inline offset(::$name) = $(off)
    end
end