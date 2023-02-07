
struct PrecessionNutationComponent{T, DA, DB}
    A::SVector{DA, T}
    B::SVector{DB, T}
    Θ₁::SVector{DB, T}
    Θ₂::SVector{DB, T}
    fun::Symbol
    nuts::Bool
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
        hasnp ? nuts[1][1:length(angles)] : nothing,
        hasnp ? nuts[2][1:length(angles)] : nothing,
        fun
    )
end

function Base.show(io::IO, p::PrecessionNutationComponent{T, DA, DB}) where {T, DA, DB}
    println(io, "PrecessionNutationComponent{$T, $DA, $DB}(nutation=$(p.nuts))")
end

struct PlanetsPrecessionNutation{T}
    ra::PrecessionNutationComponent{T}
    dec::PrecessionNutationComponent{T}
    pm::PrecessionNutationComponent{T}
end

function PlanetsPrecessionNutation(NAIFId::N,  
    data::AbstractDict{N, Dict{Symbol, Union{T, Int64, Vector{T}}}}) where {N <: Integer, T}
    
    # Find nutation coefficients
    sid = "$(NAIFId)"
    nutsid = NAIFId

    if NAIFId < 1000 && NAIFId > 100
        nutsid = parse(N, sid[1])
    end

    if NAIFId < 9
        Logging.@warn "[Basic/Orient] IAU frame for point $NAIFId does not exist - IGNORED."
    else 

        # Get nutation-precession angles, if present
        if haskey(data, nutsid) && haskey(data[nutsid], :nut_prec_angles)
            nut = data[nutsid][:nut_prec_angles] .* T(pi/180)
            nuts = (@views(nut[1:2:end]), @views(nut[2:2:end]))
        else
            nuts = nothing 
        end

        hasnuts = !isnothing(nuts)

        # Build planet precession/nutation
        return PlanetsPrecessionNutation(
            # Right ascension
            PrecessionNutationComponent{T}(hasnuts, data[NAIFId], nuts, :ra, :sin),
            # Declination 
            PrecessionNutationComponent{T}(hasnuts, data[NAIFId], nuts, :dec, :cos),
            # Polar motion
            PrecessionNutationComponent{T}(hasnuts, data[NAIFId], nuts, :pm, :sin),
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

function _add_sine!(β, B, f, i)
    push!(
        β.args, 
        Expr(
            :call, 
            :*,
            B,
            Expr(:call, f, Expr(:call, :getindex, :Θᵢ, i))
        )
    )
    nothing
end

function _iau_angles(ts::Symbol, tp::Symbol, A::NV, B::NM, Θ₁::NM, Θ₂::NM, 
    χ::Symbol, nuts::Bool, factorp=1.0, factors=1.0) where {NV<:AbstractVector, NM<:AbstractVector}
    
    # β = ∑AᵢTⁱ + ∑ Bᵢ χ(θ₀ᵢ + θ₁ᵢ t)
    δχ = χ == :cos ? :sin : :cos
    δ²χ = χ == :cos ? :cos : :sin
    sχ = χ == :cos ? -1.0 : 1.0
    s²χ = χ == :cos ? -1.0 : -1.0

    β = Expr(:call, :+, )
    δβ = Expr(:call, :+, )
    δ²β = Expr(:call, :+, )

    # polynomial
    βp = Expr(:call, :+, )
    δβp = Expr(:call, :+, )
    δ²βp = Expr(:call, :+, )

    if length(A) > 0
        for i in eachindex(A)
            Aᵢ = A[i]
            if !(Aᵢ ≈ 0.0)
                # angle 
                push!(
                    βp.args, 
                    (i-1) == 0 
                        ? Aᵢ : Expr(:call, :*, Aᵢ, Expr(:call, :^, tp, Float64(i-1)))
                )
                
                # first derivative
                δAᵢ = (i-1) * Aᵢ
                if i > 1 && !(δAᵢ ≈ 0.0)
                    push!(
                        δβp.args, 
                        (i-2) == 0 
                            ? δAᵢ : Expr(:call, :*, δAᵢ, Expr(:call, :^, tp, Float64(i-2)))
                    )

                    # second derivative
                    δ²Aᵢ = (i-1) * (i-2) * Aᵢ
                    if i > 2 && !(δ²Aᵢ ≈ 0.0)
                        push!(
                            δ²βp.args, 
                            (i-3) == 0 
                                ? δ²Aᵢ : Expr(:call, :*, δ²Aᵢ, Expr(:call, :^, tp, Float64(i-3)))
                        )
                    end    
                end
            end
        end
    end
    δβp_ = Expr(:call, :*, δβp, 1/factorp)
    δ²βp_ = Expr(:call, :*, δ²βp, 1/factorp^2)
    if length(βp.args) != 1 push!(β.args, βp) end 
    if length(δβp.args) != 1 push!(δβ.args, δβp_) end 
    if length(δ²βp.args) != 1 push!(δ²β.args, δ²βp_) end 

    # sinusoidal 
    if nuts
        βs = Expr(:call, :+, )
        δβs = Expr(:call, :+, )
        δ²βs = Expr(:call, :+, )
        for i in eachindex(B)
            Bᵢ = B[i]
            if !(Bᵢ ≈ 0.0)
                Θᵢ = Expr(:call, :+, Θ₁[i], Expr(:call, :*, Θ₂[i], ts))
                # angle
                _add_sine!(βs, Bᵢ, χ, i)

                # first derivative
                δBᵢ = sχ * Bᵢ * Θ₂[i] 
                if !(δBᵢ ≈ 0)
                    _add_sine!(δβs, δBᵢ, δχ, i)

                    # second derivative
                    δ²Bᵢ = s²χ * Bᵢ * Θ₂[i]^2 
                    if !(δ²Bᵢ ≈ 0)
                        _add_sine!(δ²βs, δ²Bᵢ, δ²χ, i)

                    end
                end
            end
        end
        δβs_ = Expr(:call, :*, δβs, 1/factors)
        δ²βs_ = Expr(:call, :*, δ²βs, 1/factors^2)
        if length(βs.args) != 1 push!(β.args, βs) end 
        if length(δβs.args) != 1 push!(δβ.args, δβs_) end 
        if length(δ²βs.args) != 1 push!(δ²β.args, δ²βs_) end 
    end
    
    β = length(β.args) == 1 ? :(0.0) : β
    δβ = length(δβ.args) == 1 ? :(0.0) : δβ
    δ²β = length(δ²β.args) == 1 ? :(0.0) : δ²β

    return β, δβ, δ²β
end

function _iau_angles(ts::Symbol, tp::Symbol,  p::PrecessionNutationComponent, fp, fs)
    return _iau_angles(ts, tp, p.A, p.B, p.Θ₁, p.Θ₂, p.fun, p.nuts, fp, fs)
end

function _iau_angles(p)
    return (
        _iau_angles(:T, :T, getproperty(p, :ra), Tempo.CENTURY2SEC, Tempo.CENTURY2SEC), 
        _iau_angles(:T, :T, getproperty(p, :dec), Tempo.CENTURY2SEC, Tempo.CENTURY2SEC),
        _iau_angles(:T, :d, getproperty(p, :pm), Tempo.DAY2SEC, Tempo.CENTURY2SEC)
    )
end

function orient_planets_angles(p, name)
    Θᵢ = _compute_thetas(:T, p.ra.Θ₁, p.ra.Θ₂)
    angles = _iau_angles(p)

    a, ad, add = angles[1]
    del, dd, ddd = angles[2]
    w, wd, wdd = angles[3]

    bname = lowercase(eval(:(String($name))))

    fname = Symbol("orient_angles_$bname")
    fdname = Symbol("orient_d_angles_$bname")
    fddname = Symbol("orient_dd_angles_$bname")

    return eval(quote
        (   
            function ($fname)(sec::Number)
                d = sec/Tempo.DAY2SEC
                T = d/Tempo.CENTURY2DAY
                Θᵢ = $Θᵢ
                return $a, $del, $w
            end,
            function ($fdname)(sec::Number)
                d = sec/Tempo.DAY2SEC
                T = d/Tempo.CENTURY2DAY
                Θᵢ = $Θᵢ
                return $a, $del, $w, $ad, $dd, $wd
            end,
            function ($fddname)(sec::Number)
                d = sec/Tempo.DAY2SEC
                T = d/Tempo.CENTURY2DAY
                Θᵢ = $Θᵢ
                return $a, $del, $w, $ad, $dd, $wd, $add, $ddd, $wdd
            end
        )
    end)
end
