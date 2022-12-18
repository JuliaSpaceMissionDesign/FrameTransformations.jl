# Points

abstract type RootPoint <: AbstractFramePoint end  
abstract type FixedPoint <: AbstractFramePoint end 
abstract type TimePoint <: AbstractFramePoint end 

abstract type EphemerisPoint <: TimePoint end
abstract type UpdatablePoint <: TimePoint end






macro point(name::Symbol, id::Int, subtype::Symbol)

    name_str = String(name)


end