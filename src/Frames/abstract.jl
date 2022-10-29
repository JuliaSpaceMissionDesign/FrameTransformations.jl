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
    AbstractUpdatableFrame

A type representing frames which can be updated as a function of the time.
Once a subtype is defined, a [`update!`](@ref) method associated to it shall be 
also defined. Any subtype of `AbstractUpdatableFrame` shall have also at least 
a field `e` where the leatest update epoch is saved.
"""
abstract type AbstractUpdatableFrame <: AbstractFrame end

"""
    update!(frame::F, args...) where {F<:AbstractUpdatableFrame}

Abstract method to update an [`AbstractUpdatableFrame`](@ref) frame at epoch.

!!! warning 
    This shall be implemented for each subtype. 
"""
function update!(::F, args...) where {F<:AbstractUpdatableFrame}
    throw(error("[Frames] `update!` method shall be defined for $F as it is updatable"))
end

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

function Base.show(io::IO, ::F) where {F<:AbstractFrame}
    println(io, join(split(String(Symbol(F)), r"(?=[A-Z])"), " ") )
end