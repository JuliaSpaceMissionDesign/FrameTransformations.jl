"""
    AbstractFrame

A type representing all reference frames.
"""
abstract type AbstractFrame end

"""
    AbstractInertialFrame

A type representing all inertial reference frames.
"""
abstract type AbstractInertialFrame <: AbstractFrame end

"""
    AbstractDynamicFrame

A type representing all non-inertial reference frames.
"""
abstract type AbstractDynamicFrame <: AbstractFrame end

"""
    AbstractFixedOffsetFrame

A type representing all reference frames which are defined as a fixed offset 
from another frame.
"""
abstract type AbstractFixedOffsetFrame <: AbstractFrame end

"""
    AbstractFrozenFrame

A type representing all reference frames which are defined freezing another frame
at a specified epoch. 
"""
abstract type AbstractFrozenFrame <: AbstractInertialFrame end


"""
    AbstractRotatingFrame

A type representing rotating frames.
"""
abstract type AbstractRotatingFrame <: AbstractDynamicFrame end


# abstract type AbstractDynamicalObjectFrame <: AbstractDynamicFrame end