using Logging
using StaticArrays
using OrderedCollections
using ReferenceFrameRotations
using BenchmarkTools 

# using Basic.Utils: TPC, load
# data = load(TPC("temp/pck00010.tpc"));
# ids = 399
# p = PlanetsPrecessionNutation(ids, data)


struct PrecessionNutationComponent{T, DA, DB}
    A::SVector{DA, T}
    B::SVector{DB, T}
    Θ₁::SVector{DB, T}
    Θ₂::SVector{DB, T}
    fun::Symbol
    nuts::Bool
end

function Base.show(io::IO, p::PrecessionNutationComponent{T, DA, DB}) where {T, DA, DB}
    println(io, "PrecessionNutationComponent{$T, $DA, $DB}(nutation=$(p.nuts))")
end

struct PlanetsPrecessionNutation{T}
    ra::PrecessionNutationComponent{T}
    dec::PrecessionNutationComponent{T}
    pm::PrecessionNutationComponent{T}
end

function PrecessionNutationComponent{T}(A, B, Θ₁, Θ₂, fun) where T
    if isnothing(B) || length(B) == 0
        B = SVector{0, T}()
        Θ₁ = SVector{0, T}()
        Θ₂ = SVector{0, T}()
    end
    return PrecessionNutationComponent{T, length(A), length(B)}(
        A, 
        B, 
        Θ₁, 
        Θ₂,
        fun,
        length(B) == 0 ? false : true
    )
end

function PrecessionNutationComponent{T}(hasnuts::Bool, dbid, nuts, prop, fun, 
    factor=T(pi/180)) where T
    
    # stick to NASA/NAIF PCK naming conventions
    pole = prop != :pm ? Symbol("pole_$prop") : :pm
    # look for precession/nutation angles
    nutprec = Symbol("nut_prec_$prop")
    hasnp = hasnuts && haskey(dbid, nutprec)

    angles = hasnp ? dbid[nutprec] .* factor : T[]

    return PrecessionNutationComponent{T}(
        dbid[pole] .* factor,
        angles,
        hasnp ? nuts[1] : nothing,
        hasnp ? nuts[2] : nothing,
        fun
    )
end

function PlanetsPrecessionNutation(bid::N, 
    data::OrderedDict{N, Dict{Symbol, Union{T, Vector{T}}}}) where {N<:Integer, T}
    
    # Find nutation coefficients
    sid = "$(bid)"
    nutsid = bid
    if bid < 1000 && bid > 100
        nutsid = parse(N, sid[1])
    end

    if bid < 9
        Logging.@warn "[Universe/Orient] Cannot orient IAU frames for $bid - IGNORED."
    else 
        # Get nutation-precession angles
        nut::Vector{T} = haskey(data[nutsid], :nut_prec_angles) ? 
            data[nutsid][:nut_prec_angles] .* T(pi/180) : nothing
        lnuts = length(nut)
        hasnuts = lnuts != 0
        nuts::Tuple{Vector{T}, Vector{T}} = hasnuts ? 
            (
                @views(nut[1 : lnuts÷2]), 
                @views(nut[lnuts÷2+1 : lnuts])
            ) : nothing

        # Build planet precession/nutation
        return PlanetsPrecessionNutation(
            # Right ascension
            PrecessionNutationComponent{T}(hasnuts, data[bid], nuts, :ra, :sin),
            # Declination 
            PrecessionNutationComponent{T}(hasnuts, data[bid], nuts, :dec, :cos),
            # Polar motion
            PrecessionNutationComponent{T}(hasnuts, data[bid], nuts, :pm, :sin),
        )
    end
    
end

function _compute_thetas(t::Symbol, Θ₁::NM, Θ₂::NM) where {NM<:AbstractVector}
    Θᵢ = Expr(:call, :SVector, )
    for i in eachindex(Θ₁)
        push!(
            Θᵢ.args, 
            Expr(:call, :+, Θ₁[i], Expr(:call, :*, Θ₂[i], t))
        )
    end
    return Θᵢ
end

function _iau_angles(t::Symbol, A::NV, B::NM, Θ₁::NM, Θ₂::NM, 
    χ::Symbol, nuts::Bool) where {NV<:AbstractVector, NM<:AbstractVector}
    
    # β = ∑AᵢTⁱ + ∑ Bᵢ χ(θ₀ᵢ + θ₁ᵢ t)
    δχ = χ == :cos ? :sin : :cos
    δ²χ = χ == :cos ? :cos : :sin
    sχ = χ == :cos ? -1.0 : 1.0
    s²χ = χ == :cos ? -1.0 : -1.0

    β = Expr(:call, :+, )
    δβ = Expr(:call, :+, )
    δ²β = Expr(:call, :+, )

    # polynomial
    if length(A) > 0
        for i in eachindex(A)
            Aᵢ = A[i]
            if !(Aᵢ ≈ 0.0)
                # angle 
                push!(
                    β.args, 
                    (i-1) == 0 
                        ? Aᵢ : Expr(:call, :*, Aᵢ, Expr(:call, :^, t, i-1))
                )
                
                # first derivative
                δAᵢ = (i-1) * Aᵢ
                if i > 1 && !(δAᵢ ≈ 0.0)
                    push!(
                        δβ.args, 
                        (i-2) == 0 
                            ? δAᵢ : Expr(:call, :*, δAᵢ, Expr(:call, :^, t, i-2))
                    )

                    # second derivative
                    δ²Aᵢ = (i-1) * (i-2) * Aᵢ
                    if i > 2 && !(δ²Aᵢ ≈ 0.0)
                        push!(
                            δ²β.args, 
                            (i-3) == 0 
                                ? δ²Aᵢ : Expr(:call, :*, δ²Aᵢ, Expr(:call, :^, t, i-3))
                        )
                    end    
                end
            end
        end
    end

    # sinusoidal 
    if nuts
        for i in eachindex(B)
            Bᵢ = B[i]
            if !(Bᵢ ≈ 0.0)
                Θᵢ = Expr(:call, :+, Θ₁[i], Expr(:call, :*, Θ₂[i], t))
                # angle
                push!(
                    β.args, 
                    Expr(
                        :call, 
                        :*,
                        Bᵢ,
                        Expr(:call, χ, Expr(:call, :getindex, :Θᵢ, i))
                    )
                )

                # first derivative
                δBᵢ = sχ * Bᵢ * Θ₂[i] 
                if !(δBᵢ ≈ 0)
                    push!(
                        δβ.args, 
                        Expr(
                            :call, 
                            :*,
                            δBᵢ,
                            Expr(:call, δχ, Expr(:call, :getindex, :Θᵢ, i))
                        )
                    )

                    # second derivative
                    δ²Bᵢ = s²χ * Bᵢ * Θ₂[i]^2 
                    if !(δ²Bᵢ ≈ 0)
                        push!(
                            δ²β.args, 
                            Expr(
                                :call, 
                                :*,
                                δ²Bᵢ,
                                Expr(:call, δ²χ, Expr(:call, :getindex, :Θᵢ, i))
                            )
                        )
                    end
                end
            end
        end
    end
    
    β = length(β.args) == 1 ? :(0.0) : β
    δβ = length(δβ.args) == 1 ? :(0.0) : δβ
    δ²β = length(δ²β.args) == 1 ? :(0.0) : δ²β

    return β, δβ, δ²β
end

function _iau_angles(t::Symbol, p::PrecessionNutationComponent)
    return _iau_angles(t, p.A, p.B, p.Θ₁, p.Θ₂, p.fun, p.nuts)
end

function iau_angles(p)
    return (
        _iau_angles(:T, getproperty(p, :ra)), _iau_angles(:T, getproperty(p, :dec)),
        _iau_angles(:T, getproperty(p, :pm))
    )
end

macro orient_iau_angles(bfname, p)

    Θᵢ = eval(:(_compute_thetas(:T, $p.ra.Θ₁, $p.ra.Θ₂)))
    angles = eval(:(iau_angles($p)))

    a = angles[1][1]
    ad = angles[1][2]
    add = angles[1][3]

    d = angles[2][1]
    dd = angles[2][2]
    ddd = angles[2][3]

    w = angles[3][1]
    wd = angles[3][2]
    wdd = angles[3][3]

    namedfun = Symbol("orient_angles_$bfname")
    dnamedfun = Symbol("orient_d_angles_$bfname")
    ddnamedfun = Symbol("orient_dd_angles_$bfname")

    nameddcm = Symbol("orient_dcm_$bfname")
    dnameddcm = Symbol("orient_ddcm_$bfname")
    ddnameddcm = Symbol("orient_dddcm_$bfname")

    return esc(quote
        function $(ddnamedfun)(T::Number)   
            Θᵢ = $Θᵢ 
            return $a, $ad, $add, $d, $dd, $ddd, $w, $wd, $wdd
        end
        function $(dnamedfun)(T::Number)   
            Θᵢ = $Θᵢ 
            return $a, $ad, $d, $dd, $w, $wd
        end
        function $(namedfun)(T::Number)   
            Θᵢ = $Θᵢ 
            return $a, $d, $w
        end
        function $(nameddcm)(T::Number)   
            ra, dec, w = $namedfun(T)
            return angle_to_dcm(π/2 + ra, π/2 - dec, w, :ZXZ)
        end
        function $(dnameddcm)(T::Number)   
            ra, rad, dec, decd, w, wd = $dnamedfun(T)
            ω = SA[rad, decd, wd]
            R = angle_to_dcm(π/2 + ra, π/2 - dec, w, :ZXZ)
            return R, DCM(ddcm(R, ω))
        end
        function $(ddnameddcm)(T::Number)   
            ra, rad, radd, dec, decd,decdd, w, wd, wdd = $ddnamedfun(T)
            # ω = SA[rad, decd, wd]
            # ωd = SA[radd, decdd, wdd]
            # R = angle_to_dcm(π/2 + ra, π/2 - dec, w, :ZXZ)
            # Rd = DCM(ddcm(R, ω))
            # FIXME: da implementare correttamente
        end
        nothing
    end)
end


