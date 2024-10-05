
# ------------------------------------------------------------------------------------------
# PROMOTIONS
# ------------------------------------------------------------------------------------------

# Returns the inner datatype of a given DCM 
dcm_eltype(::Union{DCM{T},Type{DCM{T}}}) where {T} = T

# Returns a promoted type for a given tuple of DCMs 
@generated function promote_dcm_eltype(::Union{T,Type{T}}) where {T<:Tuple}
    t = Union{}
    for i in 1:length(T.parameters)
        tmp = dcm_eltype(Base.unwrapva(T.parameters[i]))
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

"""
    Rotation{O, N}

A container to efficiently compute `O`-th order rotation matrices of type `N` between two 
set of axes. It stores the Direction Cosine Matrix (DCM) and its time derivatives up to 
the (`O`-1)-th order. Since this type is immutable, the data must be provided upon 
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

    Rotation{O}(dcms::DCM...) where O 

Create a `Rotation` object of order `O`. If the number of `dcms` is smaller than `O`, the 
remaining slots are filled with null DCMs, otherwise if the number of inputs is greater than 
`O`, only the first `O` DCMs are used. 

!!! warning 
    Usage of this constructor is not recommended as it may yield unexpected results to 
    unexperienced users. 
    
---

    Rotation{O}(u::UniformScaling{N}) where {O, N}
    Rotation{O, N}(u::UniformScaling) where {O, N}

Create an `O`-order identity `Rotation` object of type `N` with identity position rotation 
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
    Rotation{S1, N}(R::Rotation{S2}) where {S1, S2, N}

Transform a `Rotation` object of order `S2` to order `S1` and type `N`. The behaviour of 
these functions depends on the values of `S1` and `S2`: 

- `S1` < `S2`: Only the first `S1` components of `rot` are considered.
- `S1` > `S2`: The missing orders are filled with null DCMs.

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

julia> R3 = Rotation{3}(R1);

julia> R3[3]
DCM{Float64}:
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0
```

---

    Rotation(m::DCM{N}, ω::AbstractVector) where N 

Create a 2nd order `Rotation` object of type `N` to rotate between two set of axes `a` and 
`b` from a Direction Cosine Matrix (DCM) and the angular velocity vector `ω` of `b` with 
respect to `a`, expressed in `b`

"""
struct Rotation{S,T} <: AbstractArray{T,1}
    m::NTuple{S,DCM{T}}

    function Rotation(tup::NTuple{S,Any}) where {S}
        T = promote_dcm_eltype(tup)
        return new{S,T}(tup)
    end
end

""" 
    order(R::Rotation{O}) where O 

Return the rotation order O.
"""
@inline order(::Rotation{S,<:Any}) where {S} = S

# Julia API
Base.size(::Rotation{S,<:Any}) where {S} = (S,)
Base.getindex(R::Rotation, i) = R.m[i]
Base.length(::Rotation{S}) where {S} = S

# ------------------------------------------------------------------------------------------
# CONSTRUCTORS
# ------------------------------------------------------------------------------------------

# Varargs constructor
function Rotation(args::Vararg{Any,S}) where {S}
    return Rotation(args)
end

# Constructor with filter and auto-fill of missing DCMS 
@generated function Rotation{S}(dcms::DCM...) where {S}
    D = length(dcms)
    T = Expr(:call, :promote_dcm_eltype, :dcms)
    expr = Expr(:call, :tuple)
    for i in 1:min(S, D)
        push!(expr.args, Expr(:ref, :dcms, i))
    end
    for _ in 1:(S-D)
        push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, :zero, T), :I)))
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

@generated function Rotation{S1}(dcms::NTuple{S2,DCM{T}}) where {S1,S2,T}
    expr = Expr(:call, :tuple)
    for i in 1:min(S1, S2)
        push!(expr.args, Expr(:ref, :dcms, i))
    end
    for _ in 1:(S1-S2)
        push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, :zero, T), :I)))
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

# Constructor for S-order identity rotations! 
@generated function Rotation{S}(::UniformScaling{T}) where {S,T}
    expr = Expr(:call, :tuple)
    for i in 1:S
        if i == 1
            push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, :one, T), :I)))
        else
            push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, :zero, T), :I)))
        end
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

@generated function Rotation{S,T}(::UniformScaling) where {S,T}
    expr = Expr(:call, :tuple)
    for i in 1:S
        if i == 1
            push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, T, 1), :I)))
        else
            push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, T, 0), :I)))
        end
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

# Convert a Rotation to a different order
@generated function Rotation{S1}(rot::Rotation{S2,T}) where {S1,S2,T}
    expr = Expr(:call, :tuple)
    for i in 1:min(S1, S2)
        push!(expr.args, Expr(:ref, :rot, i))
    end
    for _ in 1:(S1-S2)
        push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, :zero, T), :I)))
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

# Convert a rotation to a different order and type
@generated function Rotation{S1,T}(rot::Rotation{S2}) where {S1,S2,T}

    expr = Expr(:call, :tuple)
    for i in 1:min(S1, S2)
        push!(expr.args, Expr(:., T, Expr(:tuple, Expr(:ref, :rot, i))))
    end
    for _ in 1:(S1-S2)
        push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, :zero, T), :I)))
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

function Rotation(m::DCM{T}, ω::AbstractVector) where {T}
    dm = DCM(ddcm(m, SVector(ω)))
    return Rotation((m, dm))
end

@inline Rotation{S}(rot::Rotation{S}) where {S} = rot

# ------------------------------------------------------------------------------------------
# TYPE CONVERSIONS
# ------------------------------------------------------------------------------------------

# Convert a Rotation to a tuple 
@generated function Base.Tuple(rot::Rotation{S,T}) where {S,T}

    expr = Expr(:call, :tuple)
    for j in 1:(3S)
        Oⱼ = (j - 1) ÷ 3 + 1
        for i in 1:(3S)
            Oᵢ = (i - 1) ÷ 3 + 1
            if Oⱼ > Oᵢ
                push!(expr.args, Expr(:call, :zero, T))
            else
                row = i - 3 * (Oᵢ - 1)
                col = j - 3 * (Oⱼ - 1)
                rom = Oᵢ - Oⱼ + 1
                push!(
                    expr.args,
                    Expr(
                        :ref,
                        Expr(:ref, :rot, rom), row, col
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

# ------------------------------------------------------------------------------------------
# OPERATIONS
# ------------------------------------------------------------------------------------------

# Inverse
""" 
    inv(rot::Rotation)

Compute the invese of the rotation object `rot`. The operation is efficiently performed by 
taking the transpose of each rotation matrix within `rot`.
"""
Base.inv(rot::Rotation) = _inverse_rotation(rot)

@generated function _inverse_rotation(rot::Rotation{S,T}) where {S,T}
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

# Product between Rotations
@inline Base.:*(r1::Rotation{S,<:Any}, r2::Rotation{S,<:Any}) where {S} = _compose_rotation(r1, r2)
function Base.:*(::Rotation{S1,<:Any}, ::Rotation{S2,<:Any}) where {S1,S2}
    throw(DimensionMismatch("Cannot multiply two `Rotation` types of order $S1 and $S2"))
end

@generated function _compose_rotation(A::Rotation{S,<:Any}, B::Rotation{S,<:Any}) where {S}
    expr = Expr(:call, :Rotation)
    for i in 1:S
        sum_expr = Expr(:call, :+)
        for j in 1:i
            c = binomial(i - 1, j - 1)
            aᵢ = Expr(:ref, :A, i - j + 1)
            bᵢ = Expr(:ref, :B, j)
            push!(sum_expr.args, Expr(:call, :*, c, aᵢ, bᵢ))
        end
        push!(expr.args, sum_expr)
    end

    return quote
        @inbounds $(expr)
    end
end

@inline Base.:*(A::Rotation{S}, v::Translation{S}) where {S} = _apply_rotation(A, v)
@inline Base.:*(A::Rotation{S}, v::AbstractVector{T}) where {S,T} = _apply_rotation(A, v)
@inline function Base.:*(::Rotation{S1}, ::Translation{S2}) where {S1,S2}
    throw(DimensionMismatch("Cannot apply Rotation of order $S1 to Translation of order $S2"))
end

# Product between Rotation and a Translation
@generated function _apply_rotation(R::Rotation{S,Nr}, v::Translation{S,Nv}) where {S,Nr,Nv}
    # Apply rotation on a translation vector with the same size
    #
    #                    n
    # z(n) = R(n)*v(1) + ∑ binomial(n, k) * R(n-k) * v(k)
    #                   k=1

    expr = Expr(:call, :tuple)
    for i in 1:S
        sumexpr = Expr(:call, :+)
        push!(expr.args, sumexpr)
        push!(sumexpr.args, Expr(:call, :*, Expr(:ref, :R, i), Expr(:ref, :v, 1)))

        for j in 1:i-1
            push!(
                sumexpr.args,
                Expr(
                    :call, :*,
                    binomial(i - 1, j),
                    Expr(:ref, :R, i - j),
                    Expr(:ref, :v, j + 1),
                )
            )
        end
    end
    return quote
        Base.@_inline_meta
        @inbounds Translation($(expr))
    end
end

@inline Base.:*(A::Rotation, v::SVector) = _apply_rotation(A, v)

# Compute product between Rotation and a "proper" SVector (returning a Translation)
@generated function _apply_rotation(R::Rotation{Sr,Nr}, v::SVector{Sv,Nv}) where {Sr,Sv,Nr,Nv}

    if Sv != 3Sr
        throw(
            DimensionMismatch(
                "Cannot apply Rotation of order $Sr to a $(Sv) vector",
            )
        )
    end

    expr = Expr(
        :call,
        Expr(:curly, :SVector, Sv, Nv),
        Expr(
            :call,
            :_apply_rotation, :R,
            Expr(:call, Expr(:curly, :Translation, Sr), :v)
        )
    )

    return quote
        Base.@_inline_meta
        @inbounds $expr
    end
end

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
                        :(A[$mi][$k, $w] * b[3*($j-1)+$w]) for w in 1:3
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

        T = Base.promote_op(LinearAlgebra.matprod, Na, Nb)
        $retexpr
    end
end