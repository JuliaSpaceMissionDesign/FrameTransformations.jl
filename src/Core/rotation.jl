# ----
# PROMOTIONS

# Returns the inner datatype of a given DCM 
_dcm_type(::Union{DCM{T}, Type{DCM{T}}}) where {T} = T

# Returns a promoted type for a given tuple of DCMs 
@generated function promote_dcm_eltype(::Union{T, Type{T}}) where {T<:Tuple}
    t = Union{}
    for i in 1:length(T.parameters)
        tmp = _dcm_type(Base.unwrapva(T.parameters[i]))
        t = :(promote_type($t, $tmp))
    end
    return quote
        Base.@_inline_meta
        $t
    end
end

# ----
# ROTATION TYPE DEFINITION

struct Rotation{O, N}
    m::NTuple{O, DCM{N}}

    function Rotation(x::NTuple{O, Any}) where {O}
        T = promote_dcm_eltype(x)
        return new{O, T}(x)
    end
end

@inline order(::Rotation{S,<:Any}) where {S} = S

# Julia API
Base.size(::Rotation{S,<:Any}) where {S} = (S,)
Base.getindex(R::Rotation, i) = R.m[i]
Base.length(::Rotation{S}) where S = S

# Static Arrays API 
StaticArrays.Size(::Rotation{S,N}) where {S,N} = Size((S,))
StaticArrays.similar_type(::Rotation{S,N}) where {S,N} = Rotation{S,N}
StaticArrays.similar_type(::Rotation{S,<:Any}, ::Type{N}) where {S,N} = Rotation{S,N}

# ----
# Constructors

# Generic rotation constructor 
@generated function Rotation(dcms::DCM...)
    S = length(dcms)

    expr = Expr(:call, :tuple)
    for i in 1:S 
        push!(expr.args, Expr(:ref, :dcms, i))
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

# Constructor with filter and auto-fill of missing DCMS 
@generated function Rotation{S}(dcms::DCM...) where {S}
    D = length(dcms)
    expr = Expr(:call, :tuple)
    for i in 1:min(S, D)
        push!(expr.args, Expr(:ref, :dcms, i))
    end
    for _ in 1:(S - D)
        push!(expr.args, Expr(:call, :DCM, Expr(:call, :*, 0, :I)))
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

@generated function Rotation{S1}(dcms::NTuple{S2, DCM{N}}) where {S1, S2, N}
    expr = Expr(:call, :tuple)
    for i in 1:min(S1, S2)
        push!(expr.args, Expr(:ref, :dcms, i))
    end
    for _ in 1:(S1 - S2)
        push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, N, 0), :I)))
    end
    return quote 
        @inbounds Rotation($(expr))
    end
end

# Constructor for S-order identity rotations! 
@generated function Rotation{S}(::UniformScaling{N}) where {S,N}
    expr = Expr(:call, :tuple)
    for i in 1:S 
        if i == 1
            push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, N, 1), :I)))
        else
            push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, N, 0), :I)))
        end
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

@generated function Rotation{S,N}(::UniformScaling) where {S,N}
    expr = Expr(:call, :tuple)
    for i in 1:S 
        if i == 1
            push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, N, 1), :I)))
        else
            push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, N, 0), :I)))
        end
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

# Convert a Rotation to a different order
@generated function Rotation{S1}(rot::Rotation{S2,N}) where {S1,S2,N}
    expr = Expr(:call, :tuple)
    for i in 1:min(S1, S2)
        push!(expr.args, Expr(:ref, :rot, i))
    end
    for _ in 1:(S1 - S2)
        push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, N, 0), :I)))
    end
    return quote 
        @inbounds Rotation($(expr))
    end
end

# Convert a rotation to a different order and type
@generated function Rotation{S1, N}(rot::Rotation{S2}) where {S1,S2,N}

    expr = Expr(:call, :tuple)
    for i in 1:min(S1, S2)
        push!(expr.args, Expr(:., N, Expr(:tuple, Expr(:ref, :rot, i))))
    end
    for _ in 1:(S1 - S2)
        push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, N, 0), :I)))
    end
    return quote 
        @inbounds Rotation($(expr))
    end
end

function Rotation(m::DCM{N}, ω::AbstractVector) where {N}
    dm = DCM(ddcm(m, SVector(ω)))
    return Rotation((m, dm))
end

# ---
# TYPE CONVERSIONS

# Convert a Rotation to a tuple 
@generated function Base.Tuple(rot::Rotation{S,N}) where {S,N}

    expr = Expr(:call, :tuple)
    for j in 1:(3S)
        Oⱼ = (j - 1) ÷ 3 + 1
        for i in 1:(3S)
            if Oⱼ > Oᵢ
                push!(expr.args, Expr(:call, N, 0))
            else 
                row = i - 3 * (Oᵢ - 1)
                col = j - 3 * (Oⱼ - 1)
                rom = Oᵢ - Oⱼ + 1
                push!(
                    expr.args, 
                    Expr(
                        :ref, 
                        Expr(:ref, :rot, rom), row,  col
                    )
                )
            end
        end
    end

    return quote
        Base.@_inline_meta
        @inbounds $(expr)
    end
end

# Generic Rotation-to-StaticArrays conversions
@inline function (::Type{SA})(rot::Rotation) where {SA<:StaticArray}
    return SA(Tuple(rot))
end

@inline MMatrix(rot::Rotation{S,N}) where {S,N} = MMatrix{3 * S,3 * S}(rot)
@inline SMatrix(rot::Rotation{S,N}) where {S,N} = SMatrix{3 * S,3 * S}(rot)

# ----
# OPERATIONS 

# Inverse
Base.inv(rot::Rotation) = _inverse_rotation(rot)
@generated function _inverse_rotation(rot::Rotation{S,N}) where {S,N}
    expr = Expr(:call, :Rotation,)
    for i in 1:S 
        push!(
            expr.args, Expr(:call, :adjoint, Expr(:ref, :rot, i))
        )
    end
    return quote 
        @inbounds $(expr)
    end
end

# Rotations multiplication 
@inline Base.:*(A::Rotation{S,<:Any}, B::Rotation{S,<:Any}) where {S} = _compose_rotation(A, B)
function Base.:*(::Rotation{S1,<:Any}, ::Rotation{S2,<:Any}) where {S1, S2}
    throw(DimensionMismatch("Cannot multiply two `Rotation` types of order $S1 and $S2"))
end

@generated function _compose_rotation(A::Rotation{S,<:Any}, B::Rotation{S,<:Any}) where {S}
    expr = Expr(:call, :Rotation)

    for i in 1:S
        sum_expr = Expr(:call, :+)
        for j in 1:i
            c = binomial(i - 1, j - 1)
            ai = Expr(:ref, Expr(:., :A, QuoteNode(:m)), i - j + 1)
            bi = Expr(:ref, Expr(:., :B, QuoteNode(:m)), j)

            push!(sum_expr.args, Expr(:call, :*, c, ai, bi))
        end
        push!(expr.args, sum_expr)
    end

    return quote
        @inbounds $(expr)
    end
end

@inline Base.:*(A::Rotation, v::AbstractVector) = _apply_rotation(A, v)

# Function to compute product between Rotation and a generic vector
@generated function _apply_rotation(A::Rotation{S,Na}, b::AbstractVector{Nb}) where {S,Na,Nb}
    exprs = [[Expr(:call, :+) for _ in 1:3] for _ in 1:S]

    for i in 1:S
        for j in 1:i
            for k in 1:3
                mi = i - j + 1
                push!(
                    exprs[i][k].args,
                    StaticArrays.combine_products([
                        :(A[$mi][$k, $w] * b[3 * ($j - 1) + $w]) for w in 1:3
                    ]),
                )
            end
        end
    end

    sa = 3 * S
    retexpr = :(@inbounds return similar_type(b, T, Size($sa))(tuple($((exprs...)...))))
    return quote
        Base.@_inline_meta
        length(b) != $sa && throw(
            DimensionMismatch(
                "Cannot multiply `Rotation` of size ($($sa), $($sa)) and a $(size(b)) vector",
            ),
        )

        T = Base.promote_op(matprod, Na, Nb)
        $retexpr
    end
end

# Function to compute product between Rotation and a static vector
@generated function _apply_rotation(A::Rotation{S,Na}, b::StaticVector{sb, Nb}) where {S,sb,Na,Nb}
    sa = 3 * S
    sa != sb && throw(
        DimensionMismatch(
            "Cannot multiply `Rotation` of size ($sa, $sa) and a ($sb,) length vector"
        ),
    )

    exprs = [[Expr(:call, :+) for _ in 1:3] for _ in 1:S]
    for i in 1:S
        for j in 1:i
            for k in 1:3
                push!(
                    exprs[i][k].args,
                    StaticArrays.combine_products([
                        :(A[$i - $j + 1][$k, $w] * b[3 * ($j - 1) + $w]) for w in 1:3
                    ]),
                )
            end
        end
    end

    retexpr = :(@inbounds return similar_type(b, T, Size($sa))(tuple($((exprs...)...))))
    return quote
        Base.@_inline_meta
        T = Base.promote_op(matprod, Na, Nb)
        $retexpr
    end
end
