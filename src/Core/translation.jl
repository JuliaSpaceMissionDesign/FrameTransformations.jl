
# ------------------------------------------------------------------------------------------
# PROMOTIONS
# ------------------------------------------------------------------------------------------

const SVector3{T} = SVector{3,T}

svector3_eltype(::Union{SVector3{T},Type{SVector3{T}}}) where {T} = T

# Returns a promoted type for a given tuple of SVector3
@generated function promote_svector3_eltype(::Union{T,Type{T}}) where {T<:Tuple}
    t = Union{}
    for i in 1:length(T.parameters)
        tmp = svector3_eltype(Base.unwrapva(T.parameters[i]))
        t = :(promote_type($t, $tmp))
    end
    return quote
        Base.@_inline_meta
        $t
    end
end

# ------------------------------------------------------------------------------------------
# TYPE DEF
# ------------------------------------------------------------------------------------------

struct Translation{S,T<:Number} <: AbstractArray{T,1}
    v::NTuple{S,SVector3{T}}

    function Translation(tup::NTuple{S,Any}) where {S}
        T = promote_svector3_eltype(tup)
        return new{S,T}(tup)
    end
end

function Base.summary(io::IO, ::Translation{S,N}) where {S,N}
    return print(io, "Translation{$S, $N}")
end

order(::Translation{S,<:Any}) where {S} = S

# Julia API
Base.size(::Translation{S,<:Any}) where {S} = (S,)
Base.getindex(t::Translation, i) = t.v[i]
Base.length(::Translation{S}) where {S} = S
Base.keys(t::Translation) = keys(t.v)

# ------------------------------------------------------------------------------------------
# CONSTRUCTORS
# ------------------------------------------------------------------------------------------

# Varargs constructor
function Translation(args::Vararg{<:Number,S}) where {S}
    O, r = divrem(S, 3)
    if r != 0
        throw(
            DimensionMismatch(
                "Cannot initialize a Translation with $S arguments: shall be divisible by 3.")
        )
    end
    T = promote_type(typeof.(args)...)
    return Translation(ntuple(i -> SVector3{T}(@views(args[((i-1)*3+1):(3*i)])), O))
end

# Empty constructor
@generated function Translation{S,T}() where {S,T}
    expr = Expr(:call, :Translation)
    for _ in 1:3*S
        push!(expr.args, zero(T))
    end
    return quote
        Base.@_inline_meta
        $expr
    end
end

# Empty constructor
@generated function Translation{O,T}(args::Vararg{<:Number,L}) where {O,L,T}
    expr = Expr(:call, :Translation)
    for i in 1:L
        push!(
            expr.args, Expr(:ref, :args, i)
        )
    end
    for _ in 1:(3*O-L)
        push!(expr.args, zero(T))
    end
    return quote
        Base.@_inline_meta
        $expr
    end
end


# Convert to a different order auto-fill of missing SVector3
@generated function Translation{S1}(tr::Translation{S2,T}) where {S1,S2,T}
    expr = Expr(:call, :Translation)
    for i in 1:min(S1, S2)
        v = Expr(:ref, :tr, i)
        for j in 1:3
            push!(expr.args, Expr(:ref, v, j))
        end
    end
    for _ in 1:3*(S1-S2)
        push!(expr.args, zero(T))
    end
    return quote
        Base.@_inline_meta
        @inbounds $(expr)
    end
end

@inline Translation{S}(tr::Translation{S}) where {S} = tr

# ------------------------------------------------------------------------------------------
# TYPE CONVERSIONS
# ------------------------------------------------------------------------------------------

# Convert to Tuple 
@generated function Base.Tuple(tr::Translation{S,N}) where {S,N}
    expr = Expr(:call, :tuple)
    for i in 1:S
        v = Expr(:ref, :tr, i)
        for j in 1:3
            push!(expr.args, Expr(:ref, v, j))
        end
    end
    return quote
        Base.@_inline_meta
        @inbounds $expr
    end
end

# Generic convert to SVector
@inline function (::Type{SA})(tr::Translation) where {SA<:StaticArray}
    return SA(Tuple(tr))
end

# Constructor from SVector 
@generated function Translation{St}(v::SVector{Sv,T}) where {St,Sv,T}
    if rem(Sv, 3) != 0
        throw(
            DimensionMismatch(
                "Cannot create Translation from vector with size $(Sv), shall be divisible by 3."
            )
        )
    end
    expr = Expr(:call, :Translation)
    for i in 1:min(Sv, 3 * St)
        push!(expr.args, Expr(:ref, :v, i))
    end
    for _ in 1:(3*St-Sv)
        push!(expr.args, zero(T))
    end
    return quote
        Base.@_inline_meta
        @inbounds $(expr)
    end
end

# ------------------------------------------------------------------------------------------
# OPERATIONS
# ------------------------------------------------------------------------------------------

function Base.:(==)(t1::Translation{S1}, t2::Translation{S2}) where {S1,S2}
    throw(
        DimensionMismatch("Cannot compare two Translations with different orders.")
    )
end

function Base.:(==)(t1::Translation{S}, t2::Translation{S}) where {S}
    if length(t1) != length(t2)
        return false
    end
    for i in eachindex(t1)
        if t1[i] != t2[i]
            return false
        end
    end
    return true
end

@generated function Base.:(+)(t1::Translation{S1}, t2::Translation{S2}) where {S1,S2}
    expr = Expr(:call, :tuple)
    if S1 ≥ S2
        for i in 1:S2
            push!(
                expr.args, Expr(:call, :(+), Expr(:ref, :t1, i), Expr(:ref, :t2, i))
            )
        end
        for i in S2+1:S1
            push!(expr.args, Expr(:ref, :t1, i))
        end
    else
        for i in 1:S1
            push!(
                expr.args, Expr(:call, :(+), Expr(:ref, :t1, i), Expr(:ref, :t2, i))
            )
        end
        for i in S1+1:S2
            push!(expr.args, Expr(:ref, :t2, i))
        end
    end
    trexpr = Expr(:call, :Translation, expr)
    return quote
        Base.@_inline_meta
        @inbounds $trexpr
    end

end

@generated function Base.:(-)(t::Translation{S}) where {S}
    expr = Expr(:call, :tuple)
    for i in 1:S
        push!(expr.args, Expr(:call, :(-), Expr(:ref, :t, i)))
    end

    trexpr = Expr(:call, :Translation, expr)
    return quote
        Base.@_inline_meta
        @inbounds $trexpr
    end
end

@generated function Base.:(-)(t1::Translation{S1}, t2::Translation{S2}) where {S1,S2}
    expr = Expr(:call, :tuple)
    if S1 ≥ S2
        for i in 1:S2
            push!(
                expr.args, Expr(:call, :(-), Expr(:ref, :t1, i), Expr(:ref, :t2, i))
            )
        end
        for i in S2+1:S1
            push!(expr.args, Expr(:ref, :t1, i))
        end
    else
        for i in 1:S1
            push!(
                expr.args, Expr(:call, :(-), Expr(:ref, :t1, i), Expr(:ref, :t2, i))
            )
        end
        for i in S1+1:S2
            push!(expr.args, Expr(:ref, :t2, i))
        end
    end
    trexpr = Expr(:call, :Translation, expr)
    return quote
        Base.@_inline_meta
        @inbounds $trexpr
    end

end
