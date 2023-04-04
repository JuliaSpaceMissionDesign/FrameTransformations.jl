"""
    AbstractEpochOrigin

Abstract type for all epoch origins.
"""
abstract type AbstractEpochOrigin end

const EPOCH_ORIGIN = (
    :JulianDate, :ModifiedJulianDate, :JulianDate2000, :ModifiedJulianDate2000
)

const EPOCH_ORIGIN_ACRONYMS = (:JD, :MJD, :J2000, :MJ2000)

const EPOCH_STARTS = (
    "-4712-01-01T12:00", "1858-11-17T00:00", "2000-01-01T12:00", "2000-01-01T00:00"
)

const EPOCH_ORIGIN_TO_J2000 = (DJ2000, DMJD, 0.0, -0.5)

for (name, acr, off, start) in
    zip(EPOCH_ORIGIN, EPOCH_ORIGIN_ACRONYMS, EPOCH_ORIGIN_TO_J2000, EPOCH_STARTS)
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

        Offset in days to shift J2000 epochs (with origin at `2000-01-01T12:00`) 
        to [`$($name_str)`](@ref) (with origin at `$($start)`)
        """
        @inline offset(::$name) = $(off)
    end
end
