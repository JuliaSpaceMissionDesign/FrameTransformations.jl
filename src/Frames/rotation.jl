export Rotation, order

import LinearAlgebra: matprod, UniformScaling
import StaticArrays: similar_type, Size, MMatrix, SMatrix

# -------------------------------------
# TYPES
# -------------------------------------

"""
    Rotation{S, N}

A container to efficiently compute `S`-th order rotation matrices of type `N` between two 
set of axes. It stores the Direction Cosine Matrix (DCM) and its time derivatives up to 
the (`S`-1)-th order. Since this type is immutable, the data must be provided upon 
construction and cannot be mutated later.

The rotation of state vector between two set of axes is computed with an ad-hoc overload 
of the product operator. For example, a 3rd order Rotation object `R`, constructed from the 
DCM `A` and its time derivatives `δA` and `δ²A` rotates a vector `v` = `[p, v, a]` as: 

`̂v = [A*p, δA*p + A*v, δ²A*p + 2δA*v + A*a]`

A `Rotation` object `R` call always be converted to a `SMatrix` or a `MMatrix` by invoking 
the proper constructor. 

### Examples 
```julia-repl 
julia> A = angle_to_dcm(π/3, :Z)
DCM{Float64}:
  0.5       0.866025  0.0
 -0.866025  0.5       0.0
  0.0       0.0       1.0

julia> R = Rotation(A);

julia> SM = SMatrix(R)
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
  0.5       0.866025  0.0
 -0.866025  0.5       0.0
  0.0       0.0       1.0

julia> MM = MMatrix(R)
3×3 MMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
  0.5       0.866025  0.0
 -0.866025  0.5       0.0
  0.0       0.0       1.0
```

---

    Rotation(dcms::DCM...)

Create a `Rotation` object from a Direction Cosine Matrix (DCM) and any of its time  
derivatives. The rotation order is inferred from the number of inputs, while the rotation 
type is obtained by promoting the DCMs types.

### Examples 
```julia-repl
julia> A = angle_to_dcm(π/3, :Z); 

julia> δA = DCM(0.0I); 

julia> δ²A = DCM(0.0I); 

julia> R = Rotation(A, δA, δ²A); 

julia> typeof(R) 
Rotation{3, Float64}

julia> R2 = Rotation(A, B, C, DCM(0.0im*I)); 

julia> typeof(R2)
Rotation{4, ComplexF64}
```

---

    Rotation{S}(u::UniformScaling{N}) where {S, N}
    Rotation{S, N}(u::UniformScaling) where {S, N}

Create an `S`-order identity `Rotation` object of type `N` with identity position rotation 
and null time derivatives.

### Examples 
```julia-repl
julia> Rotation{1}(1.0I) 
Rotation{1, Float64}(([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],))

julia> Rotation{1, Int64}(I)
Rotation{1, Int64}(([1 0 0; 0 1 0; 0 0 1],))
```

---

    Rotation{S1}(rot::Rotation{S2, N}) where {S1, S2, N}

Transform a `Rotation` object of order `S2` to order `S1`. This conversion is only possible 
if `S1` < `S2`.

### Examples 
```julia-repl
julia> A = angle_to_dcm(π/3, :Z);

julia> B = angle_to_dcm(π/4, π/6, :XY);

julia> R1 = Rotation(A, B);

julia> order(R1)
2

julia> R2 = Rotation{1}(R1);

julia> order(R2)
1

julia> R2[1] == A 
true

julia> R3 = Rotation{3}(R1)
ERROR: DimensionMismatch: Cannot convert a `Rotation` of order 2 to order 3
[...]
```

---

    Rotation(m::DCM{N}, ω::AbstractVector) where N 

Create a 2nd order `Rotation` object of type `N` to rotate between two set of axes `a` and 
`b` from a Direction Cosine Matrix (DCM) and the angular velocity vector `ω` of `b` with 
respect to `a`, expressed in `b`

---

    Rotation{S}(dcms::DCM...) where S 

Create a `Rotation` object of order `S`. If the number of `dcms` is smaller than `S`, the 
remaining slots are filled with null DCMs, otherwise if the number of inputs is greater than 
`S`, only the first `S` DCMs are used. 

!!! warning 
    Usage of this constructor is not recommended as it may yield unexpected results to 
    unexperienced users. 

### See also 
See also [`rotation3`](@ref), [`rotation6`](@ref) and [`rotation9`](@ref).
"""
struct Rotation{S,N}
    m::NTuple{S,DCM{N}}

    function Rotation(x::NTuple{S,Any}) where {S}
        T = promote_dcm_eltype(x)
        return new{S,T}(x)
    end
end

""" 
    order(R::Rotation{S}) where S 

Return the rotation order S.
"""
order(::Rotation{S,<:Any}) where {S} = S

# Julia APIs 
Base.size(::Rotation{S,<:Any}) where {S} = (3S, 3S)
Base.getindex(R::Rotation, i) = R.m[i]

function Base.summary(io::IO, ::Rotation{S,N}) where {S,N}
    return print(io, "Rotation{$S, $N}")
end

# ---
# Constructors 

# Generic rotation constructor 
@generated function Rotation(dcms::DCM...)
    S = length(dcms)
    expr = :(tuple($([Expr(:ref, :dcms, i) for i in 1:S]...)))
    return quote
        @inbounds Rotation($(expr))
    end
end

# Constructor with filter and auto-fill of missing DCMS 
@generated function Rotation{S}(dcms::DCM...) where {S}
    nd = length(dcms)
    expr = :(tuple(
        $([Expr(:ref, :dcms, i) for i in 1:min(nd, S)]...),
        $([DCM(0I) for i in 1:(S - nd)]...),
    ))
    return quote
        @inbounds Rotation($(expr))
    end
end

# Constructor for S-order identity rotations! 
@generated function Rotation{S}(::UniformScaling{N}) where {S,N}
    expr = :(tuple($([i == 1 ? DCM(N(1)I) : DCM(N(0)I) for i in 1:S]...)))
    return quote
        @inbounds Rotation($(expr))
    end
end

@generated function Rotation{S,N}(::UniformScaling) where {S,N}
    expr = :(tuple($([i == 1 ? DCM(N(1)I) : DCM(N(0)I) for i in 1:S]...)))
    return quote
        @inbounds Rotation($(expr))
    end
end

# Convert a Rotation to one with a smaller order! 
function Rotation{S1}(rot::Rotation{S2,N}) where {S1,S2,N}
    S1 > S2 &&
        throw(DimensionMismatch("Cannot convert a `Rotation` of order $S2 to order $S1"))
    return Rotation(rot.m[1:S1])
end

function Rotation(m::DCM{N}, ω::AbstractVector) where {N}
    dm = DCM(ddcm(m, SVector(ω)))
    return Rotation((m, dm))
end

# ---
# Type Conversions and Promotions 

# Static Arrays API 
Size(::Rotation{S,N}) where {S,N} = Size((3 * S, 3 * S))

similar_type(::Rotation{S,N}) where {S,N} = Rotation{S,N}
similar_type(::Rotation{S,<:Any}, ::Type{N}) where {S,N} = Rotation{S,N}

# Convert a Rotation to a tuple 
@generated function Base.Tuple(rot::Rotation{S,N}) where {S,N}
    args = []
    for j in 1:(3S)
        Oⱼ = (j - 1) ÷ 3 + 1
        for i in 1:(3S)
            Oᵢ = (i - 1) ÷ 3 + 1
            if Oⱼ > Oᵢ
                push!(args, N(0))
            else
                row = i - 3 * (Oᵢ - 1)
                col = j - 3 * (Oⱼ - 1)

                push!(args, :(rot[$Oᵢ - $Oⱼ + 1][$row, $col]))
            end
        end
    end

    return quote
        Base.@_inline_meta
        @inbounds tuple($(args...))
    end
end

# Generic Rotation-to-StaticArrays conversions
@inline function (::Type{SA})(rot::Rotation) where {SA<:StaticArray}
    return SA(Tuple(rot))
end

@inline MMatrix(rot::Rotation{S,N}) where {S,N} = MMatrix{3 * S,3 * S}(rot)
@inline SMatrix(rot::Rotation{S,N}) where {S,N} = SMatrix{3 * S,3 * S}(rot)

# Returns the inner datatype of a given DCM 
_dcm_type(::Union{DCM{T},Type{DCM{T}}}) where {T} = T

# Returns a promoted type for a given tuple of DCMs 
@generated function promote_dcm_eltype(::Union{T,Type{T}}) where {T<:Tuple}
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

# -------------------------------------
# OPERATIONS 
# -------------------------------------

# ---
""" 
    inv(rot::Rotation)

Compute the invese of the rotation object `rot`. The operation is efficiently performed by 
taking the transpose of each rotation matrix within `rot`.
"""
Base.inv(rot::Rotation) = _inverse_rot(rot)
@generated function _inverse_rot(rot::Rotation{S,N}) where {S,N}
    quote
        @inbounds Rotation($((:(rot.m[$i]') for i in 1:S)...))
    end
end

# ---
# Compose rotations 
@inline Base.:*(A::Rotation{S,<:Any}, B::Rotation{S,<:Any}) where {S} = _compose_rot(A, B)

function Base.:*(::Rotation{S1,<:Any}, ::Rotation{S2,<:Any}) where {S1,S2}
    throw(DimensionMismatch("Cannot multiply two `Rotation` types of order $S1 and $S2"))
end

@generated function _compose_rot(A::Rotation{S,<:Any}, B::Rotation{S,<:Any}) where {S}
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

# ---
# Apply rotation

@inline Base.:*(A::Rotation, b::AbstractVector) = _apply_rot(A, b)

# Function to compute product between Rotation and a generic vector
@generated function _apply_rot(A::Rotation{S,Na}, b::AbstractVector{Nb}) where {S,Na,Nb}
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
@generated function _apply_rot(A::Rotation{S,Na}, b::StaticVector{sb,Nb}) where {S,sb,Na,Nb}
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
