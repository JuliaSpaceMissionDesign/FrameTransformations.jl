module ChainRulesCoreExt

using FrameTransformations
using ChainRulesCore
using StaticArrays: SVector, SMatrix

using FrameTransformations: Rotation
using ReferenceFrameRotations: DCM
using JSMDUtils.Autodiff: derivative

# ==========================================================================================
# Helper: extract numeric values from cotangents
# ==========================================================================================

# Cotangents may arrive as plain arrays (Zygote) or as ChainRulesCore.Tangent structs
# (Mooncake via @from_chainrules). We need to handle both transparently.

# For vectors: convert cotangent Δ to a concrete SVector for dot product computation.
_to_svec(x::SVector{N}) where {N} = x
_to_svec(x::AbstractVector) = SVector(x...)
function _to_svec(x::ChainRulesCore.Tangent{<:SVector{N}}) where {N}
    return SVector(x.data...)
end
function _to_svec(x::ChainRulesCore.Tangent{<:Any, <:NamedTuple{(:data,)}})
    d = x.data
    if d isa ChainRulesCore.Tangent
        return SVector(d.backing...)
    elseif d isa Tuple
        return SVector(d...)
    else
        return SVector(d...)
    end
end
_to_svec(x) = SVector(collect(x)...)

# Scalar dot product between derivative J and cotangent Δ
function _vec_pullback_dt(J, Δ)
    Δv = _to_svec(unthunk(Δ))
    return sum(J .* Δv)
end

# ==========================================================================================
# Helper: Rotation pullback — Frobenius inner product between J and Δ
# ==========================================================================================

# derivative() on rotation3/6/9/12 returns a Vector{DCM{T}} (S elements for order S).
# The cotangent Δ may be a Rotation, a Tangent{Rotation}, or similar.
# We compute ∂L/∂t = Σᵢ ⟨Jᵢ, Δᵢ⟩_F  (Frobenius inner product per DCM)

# Extract the flat NTuple{9} data from a DCM or its tangent representation
_dcm_flat(x::DCM) = x.data  # NTuple{9,T}
_dcm_flat(x::SMatrix{3,3}) = Tuple(x)
_dcm_flat(x::AbstractMatrix) = Tuple(x)
# Mooncake tangent for DCM: Tangent{DCM, @NamedTuple{data::NTuple{9,T}}}
function _dcm_flat(x::ChainRulesCore.Tangent)
    d = x.data
    if d isa Tuple
        return d
    elseif d isa ChainRulesCore.Tangent
        return Tuple(d.backing)
    else
        return Tuple(d)
    end
end

# Frobenius inner product between two flat tuples
_frob_dot(a::NTuple{N}, b::NTuple{N}) where {N} = sum(a .* b)
_frob_dot(a::NTuple{N}, b) where {N} = sum(a .* Tuple(b))
_frob_dot(a, b) = sum(Tuple(a) .* Tuple(b))

# Extract the tuple of DCMs from the rotation cotangent
_get_rot_dcm_tuple(R::Rotation) = R.m
function _get_rot_dcm_tuple(Δ::ChainRulesCore.Tangent)
    m = Δ.m
    if m isa Tuple
        return m
    elseif m isa ChainRulesCore.Tangent
        return Tuple(m.backing)
    else
        return Tuple(m)
    end
end
_get_rot_dcm_tuple(Δ) = Δ.m

function _rotation_pullback_dt(J::AbstractVector, Δ)
    Δu = unthunk(Δ)
    Δm = _get_rot_dcm_tuple(Δu)
    dt = 0.0
    for i in eachindex(J)
        dt += _frob_dot(_dcm_flat(J[i]), _dcm_flat(unthunk(Δm[i])))
    end
    return dt
end

# ==========================================================================================
# rrules for vector functions: vector3, vector6, vector9, vector12
# ==========================================================================================

for (vfun, N) in (
    (:vector3, 3), (:vector6, 6), (:vector9, 9), (:vector12, 12)
)
    @eval begin
        function ChainRulesCore.rrule(::typeof($vfun), fr::FrameSystem, from, to, ax, t::Number)
            val = $vfun(fr, from, to, ax, t)
            function $(Symbol(vfun, :_pullback))(Δ)
                J = derivative(τ -> $vfun(fr, from, to, ax, τ), t)
                dt = _vec_pullback_dt(J, Δ)
                return NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), dt
            end
            return val, $(Symbol(vfun, :_pullback))
        end
    end
end

# ==========================================================================================
# rrules for rotation functions: rotation3, rotation6, rotation9, rotation12
# ==========================================================================================

for rfun in (:rotation3, :rotation6, :rotation9, :rotation12)
    @eval begin
        function ChainRulesCore.rrule(::typeof($rfun), fr::FrameSystem, from, to, t::Number)
            val = $rfun(fr, from, to, t)
            function $(Symbol(rfun, :_pullback))(Δ)
                J = derivative(τ -> $rfun(fr, from, to, τ), t)
                dt = _rotation_pullback_dt(J, Δ)
                return NoTangent(), NoTangent(), NoTangent(), NoTangent(), dt
            end
            return val, $(Symbol(rfun, :_pullback))
        end
    end
end

# ==========================================================================================
# rrules for direction functions: direction3, direction6, direction9, direction12
# ==========================================================================================

for (dfun, N) in (
    (:direction3, 3), (:direction6, 6), (:direction9, 9), (:direction12, 12)
)
    @eval begin
        function ChainRulesCore.rrule(::typeof($dfun), fr::FrameSystem, name::Symbol, ax, t::Number)
            val = $dfun(fr, name, ax, t)
            function $(Symbol(dfun, :_pullback))(Δ)
                J = derivative(τ -> $dfun(fr, name, ax, τ), t)
                dt = _vec_pullback_dt(J, Δ)
                return NoTangent(), NoTangent(), NoTangent(), NoTangent(), dt
            end
            return val, $(Symbol(dfun, :_pullback))
        end
    end
end

end # module
