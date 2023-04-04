export interpolate

"""
    AbstractInterpolationMethod

Abstract type representing any interpolation method.
"""
abstract type AbstractInterpolationMethod end

"""
    interpolate(::AbstractInterpolationMethod, x)

Abstract interpolator call method.
"""
function interpolate(::AbstractInterpolationMethod, x) end
