export add_axes_bcrtod!, 
       add_axes_bcifix!, 
       add_axes_bci2000!


# Exported functions
# ==============================================

"""
    add_axes_bcrtod!(frames, axes::AbstractFrameAxes, center, data)

Add `axes` as a set of Body-Centered Rotating (BCR), True-of-Date (TOD) axes to the 
`frames` system. The center point (i.e., the reference body) is `center`. `data` is a dictionary 
containing a parsed `TPC` file. These axes are the equivalent of SPICE's `IAU_<BODY_NAME>` frames.

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(AXESID_ICRF)). If the 
    ICRF is not defined in `frames`, an error is thrown.

----

    add_axes_bcrtod!(frames, name::Symbol, axesid::Int, cname::Symbol, cid::Int, data)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro and of an [`AbstractFramePoint`](@ref) via the [`@point`](@ref)
macro.

### See also 
See also [`add_axes_rotating!`](@ref), [`add_axes_bci2000!`](@ref) and [`AXESID_ICRF`](@ref).
"""
function add_axes_bcrtod!(
    frames::FrameSystem, axes::AbstractFrameAxes, center, data
)
    return add_axes_bcrtod!(
        frames, axes_name(axes), axes_id(axes), point_id(center), data
    )
end

function add_axes_bcrtod!(
    frames::FrameSystem, name::Symbol, axesid::Int, cid::Int, data
)
    if !(has_axes(frames, AXESID_ICRF))
        throw(
            ErrorException(
                "Body-Centered Rotating (BCR), True-of-Date (TOD) axes can only be defined" * 
                " w.r.t. the ICRF (ID = $(AXESID_ICRF)), which is not defined in" * 
                " the current frames graph."
            )
        )
    end

    p = PlanetaryRotationalElements(cid, data)
    ψ, α, δ, W = build_iau_series(p)

    ex = :(
        function (T::Number)
            d = T * 36525
            ψ = $ψ
            α = $α 
            δ = $δ 
            w = $W
            return angle_to_dcm(π / 2 + α, π / 2 - δ, w, :ZXZ)
        end
    )

    f = @RuntimeGeneratedFunction(ex)

    # TODO: implement higher order derivatives
    # Insert new axes in the frame system 
    return add_axes_rotating!(
        frames, name, axesid, AXESID_ICRF, f
    )

end


"""
    add_axes_bci2000!(frames, axes::AbstractFrameAxes, center, data)

Add `axes` as a set of Body-Centered Inertial (BCI) axes at J2000 relative to the `frames` 
system. The center point (i.e., the reference body) is `center` and can either be the point 
ID or an [`AbstractFramePoint`](@ref) instance. `data` is a dictionary containing a parsed 
`TPC` file. 

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(AXESID_ICRF)). If the 
    ICRF is not defined in `frames`, an error is thrown.

----

    add_axes_bci2000!(frames, name::Symbol, axesid::Int, cid::Int, data)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.
    
### See also 
See also [`add_axes_bcifix!`](@ref), [`add_axes_bcrtod!`](@ref) and [`AXESID_ICRF`](@ref).
"""
function add_axes_bci2000!(
    frames::FrameSystem, axes::AbstractFrameAxes, center, data
)
    return add_axes_bci2000!(
        frames, axes_name(axes), axes_id(axes), point_id(center), data
    )
end

function add_axes_bci2000!(
    frames::FrameSystem, name::Symbol,  axesid::Int,  cid::Int,  data
)
    # Insert new axes in the frame system 
    return add_axes_bcifix!(frames, name, axesid, cid, 0, data)

end


"""
    add_axes_bcifix!(frames, axes::AbstractFrameAxes, center, e::Epoch, data)

Add `axes` as a set of Body-Centered Inertial (BCI) at epoch `e`, relative to the `frames` 
system. The center point (i.e., the reference body) is `center` and can either be the point 
ID or an [`AbstractFramePoint`](@ref) instance. `data` is a dictionary containing a parsed 
`TPC` file. 

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(AXESID_ICRF)). If the 
    ICRF is not defined in `frames`, an error is thrown.

----

    add_axes_bcifix!(frames, name::Symbol, axesid::Int, cid::Int, epoch::Number, data)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.
    
### See also 
See also [`add_axes_bci2000!`](@ref), [`add_axes_bcrtod!`](@ref) and [`AXESID_ICRF`](@ref).
"""
function add_axes_bcifix!(
    frames::FrameSystem, axes::AbstractFrameAxes, center, epoch::Epoch, data
)
    return add_axes_bcifix!(
        frames, axes_name(axes), axes_id(axes), point_id(center), 
        j2000s(epoch), data
    )
end

function add_axes_bcifix!(
    frames::FrameSystem, name::Symbol, axesid::Int, cid::Int, epoch::Number, data
)   
    if !(has_axes(frames, AXESID_ICRF))
        throw(
            ErrorException(
                "Body-Centered Inertial at Epoch axes can only be defined" * 
                " w.r.t. the ICRF (ID = $(AXESID_ICRF)), which is not defined in" * 
                " the current frames graph."
            )
        )
    end

    p = PlanetaryRotationalElements(cid, data)
    ψ, α, δ, W = build_iau_series(p)

    ex = :(
        function (T::Number)
            d = T * 36525
            ψ = $ψ
            α = $α 
            δ = $δ 
            w = $W
            return angle_to_dcm(π / 2 + α, π / 2 - δ, w, :ZXZ)
        end
    )

    f = @RuntimeGeneratedFunction(ex)
    dcm = f(epoch / Tempo.CENTURY2SEC)
    
    # Insert new axes in the frame system 
    return add_axes_fixedoffset!(frames, name, axesid, AXESID_ICRF, dcm)
end



# Low-level routines
# ==============================================

struct PrecNutComponent{T}
    A::Vector{T}
    B::Vector{T}
    Θ::Vector{T}
    fun::Symbol
    nuts::Bool
    max_deg::Int
end

function PrecNutComponent{T}(A, B, Θ, fun, max_phase_deg=1) where {T}
    if isnothing(B) || length(B) == 0
        B = []
        Θ = []
    end
    return PrecNutComponent{T}(
        A, B, Θ, fun, length(B) == 0 ? false : true, max_phase_deg
    )
end

function PrecNutComponent{T}(
    hasnuts::Bool, dbid, nuts, prop, fun, max_phase_deg, factor=T(pi/180)
) where {T}
    # stick to NASA/NAIF PCK naming conventions
    pole = prop != :pm ? Symbol("pole_$prop") : :pm
    
    # look for precession/nutation angles
    nutprec = Symbol("nut_prec_$prop")
    hasnuts = hasnuts && haskey(dbid, nutprec)
    angles = hasnuts ? dbid[nutprec] .* factor : T[]

    return PrecNutComponent{T}(
        dbid[pole] .* factor, 
        angles, 
        hasnuts ? nuts : nothing, 
        fun, 
        max_phase_deg
    )
end

struct PlanetaryRotationalElements{T}
    ra::PrecNutComponent{T}
    dec::PrecNutComponent{T}
    pm::PrecNutComponent{T}
end

function PlanetaryRotationalElements(NAIFId::N, 
    data::AbstractDict{Int, Dict{Symbol, Vector{T}}}) where {N<:Integer, T}

    # Find nutation coefficients
    sid = "$(NAIFId)"
    nutsid = NAIFId

    if NAIFId < 1000 && NAIFId > 100
        nutsid = parse(N, sid[1])
    end

    if NAIFId < 9
        throw(
            ErrorException("IAU frame for point $NAIFId does not exist.")
        )
    else

        # Get nutation-precession angles, if present
        if haskey(data, nutsid) && haskey(data[nutsid], :nut_prec_angles)
            nuts = data[nutsid][:nut_prec_angles] .* T(pi/180)
        else
            nuts = nothing
        end
        hasnuts = !isnothing(nuts)

        # Get maximum phase degree
        if haskey(data, nutsid) && haskey(data[nutsid], :max_phase_degree)
            max_phase_deg = convert(Int64, data[nutsid][:max_phase_degree])
        else
            max_phase_deg = 1
        end

        # Build planet precession/nutation components
        return PlanetaryRotationalElements(
            # Right ascension
            PrecNutComponent{T}(
                hasnuts, data[NAIFId], nuts, :ra, :sin, max_phase_deg
            ),
            # Declination 
            PrecNutComponent{T}(
                hasnuts, data[NAIFId], nuts, :dec, :cos, max_phase_deg
            ),
            # Prime meridian
            PrecNutComponent{T}(
                hasnuts, data[NAIFId], nuts, :pm, :sin, max_phase_deg
            ),
        )
    end
end

function _compute_poly(t, A)
    β = Expr(:call, :+)
    @inbounds for i in eachindex(A)
        Aᵢ = A[i]
        if !(Aᵢ ≈ 0.0)
            push!(
                β.args,
                (i - 1) == 0 ? Aᵢ : Expr(:call, :*, Aᵢ, Expr(:call, :^, t, Int(i - 1))),
            )
        end
    end
    if length(β.args) == 1
        return :(0.0)
    end
    return β
end

function _add_sincos!(β, B, f, i, vec)
    push!(
        β.args, 
        Expr(
            :call, 
            :*, 
            B, 
            Expr(
                :call, 
                f, 
                Expr(:call, :getindex, vec, i)
            )
        )
    )
    return nothing
end

function _iau_sine(B, χ, angle)
    βs = Expr(:call, :+)
    for i in eachindex(B)
        Bᵢ = B[i]
        if !(Bᵢ ≈ 0.0)
            _add_sincos!(βs, Bᵢ, χ, i, angle)
        end
    end
    return βs
end

function _iau_angle(tp::Symbol, A, B, χ, nuts)
    # β = ∑ Aᵢ⋅tpⁱ + ∑ Bᵢ χ(∑ θᵢⱼ⋅tsʲ)

    # Now this is written as:
    # β = ∑ᵢ Aᵢ ⋅ tpⁱ + ∑ⱼBⱼ ⋅ χ(∑ₖψₖ)
    # where here ψ = θ⋅ts

    β = Expr(:call, :+)
    # polynomial 
    if length(A) > 0
        βp = _compute_poly(tp, A)
        if length(βp.args) > 1
            push!(β.args, βp)
        end
    end

    # sinusoidal
    if nuts && !(all(B .== 0.0))
        βs = _iau_sine(B, χ, :ψ)
        if length(βs.args) > 1
            push!(β.args, βs)
        end
    end

    if length(β.args) ≤ 1
        β = :(0.0)
    end

    return β
end

function _iau_angle(tp::Symbol, p::PrecNutComponent)
    return _iau_angle(tp, p.A, p.B, p.fun, p.nuts)
end

function _compute_ψs(t, Θm, n)
    ψ = Expr(:call, :SVector)
    for i in 1:n
        ψᵢ = _compute_poly(t, @views(Θm[i, 1:end]))
        push!(ψ.args, ψᵢ)
    end
    return ψ
end

function _compute_ψs(t, p)
    n = length(p.Θ) ÷ (p.max_deg + 1)
    Θm = reshape(p.Θ, (p.max_deg + 1, n))'
    return _compute_ψs(t, Θm, n)
end

function build_iau_series(p)
    return (
        _compute_ψs(:T, p.ra),
        _iau_angle(:T, getproperty(p, :ra)),
        _iau_angle(:T, getproperty(p, :dec)),
        _iau_angle(:d, getproperty(p, :pm)),
    )
end

