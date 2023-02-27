struct PrecessionNutationComponent{T}
    A::Vector{T}
    B::Vector{T}
    Θ::Vector{T}
    fun::Symbol 
    nuts::Bool
    max_phase_deg::Int
end

function PrecessionNutationComponent{T}(A, B, Θ, fun, max_phase_deg=1) where T
    if isnothing(B) || length(B) == 0
        B = []
        Θ = []
    end
    return PrecessionNutationComponent{T}(
        A, 
        B, 
        Θ,
        fun,
        length(B) == 0 ? false : true,
        max_phase_deg
    )
end

function PrecessionNutationComponent{T}(hasnuts::Bool, dbid, nuts, prop, fun, max_phase_deg,
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
        hasnp ? nuts : nothing,
        fun,
        max_phase_deg
    )
end

function Base.show(io::IO, p::PrecessionNutationComponent{T}) where {T}
    println(io, "PrecessionNutationComponent{$T, $(length(p.A)), $(length(p.B))}(nutation=$(p.nuts))")
end

struct PlanetsPrecessionNutation{T}
    ra::PrecessionNutationComponent{T}
    dec::PrecessionNutationComponent{T}
    pm::PrecessionNutationComponent{T}
end

function PlanetsPrecessionNutation(NAIFId::N,  
    data::AbstractDict{N, Dict{Symbol, Union{T, Int, Vector{T}}}}) where {N <: Integer, T}

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
            nuts = data[nutsid][:nut_prec_angles] .* T(pi/180)
        else
            nuts = nothing 
        end

        if haskey(data, nutsid) && haskey(data[nutsid], :max_phase_degree)
            max_phase_deg = convert(Int64, data[nutsid][:max_phase_degree])
        else
            max_phase_deg = 1 
        end

        hasnuts = !isnothing(nuts)

        # Build planet precession/nutation
        return PlanetsPrecessionNutation(
            # Right ascension
            PrecessionNutationComponent{T}(hasnuts, data[NAIFId], nuts, :ra, :sin, max_phase_deg),
            # Declination 
            PrecessionNutationComponent{T}(hasnuts, data[NAIFId], nuts, :dec, :cos, max_phase_deg),
            # Polar motion
            PrecessionNutationComponent{T}(hasnuts, data[NAIFId], nuts, :pm, :sin, max_phase_deg),
        )
    end
end

function _compute_poly(t, A)
    β = Expr(:call, :+, )
    for i in eachindex(A)
        Aᵢ = A[i]
        if !(Aᵢ ≈ 0.0)
            push!(
                β.args, 
                (i-1) == 0 
                    ? Aᵢ : Expr(:call, :*, Aᵢ, Expr(:call, :^, t, Int(i-1)))
            )
        end
    end
    if length(β.args) == 1
        return :(0.0)
    end
    return β
end

function _compute_δpoly(t, A)
    δβ = Expr(:call, :+, )
    for i in eachindex(A)
        δAᵢ = (i-1) * A[i]
        if i > 1 && !(δAᵢ ≈ 0.0)
            push!(
                δβ.args, 
                (i-2) == 0 
                    ? δAᵢ : Expr(:call, :*, δAᵢ, Expr(:call, :^, t, Int(i-2)))
            )
        end
    end
    if length(δβ.args) == 1
        return :(0.0)
    end
    return δβ 
end

function _compute_δ²poly(t, A)
    δ²β = Expr(:call, :+, )
    for i in eachindex(A)
        δ²Aᵢ = (i-1) * (i-2) * A[i]
        if i > 2 && !(δ²Aᵢ ≈ 0.0)
            push!(
                δ²β.args, 
                (i-3) == 0 
                    ? δ²Aᵢ : Expr(:call, :*, δ²Aᵢ, Expr(:call, :^, t, Int(i-3)))
            )
        end    
    end
    if length(δ²β.args) == 1
        return :(0.0)
    end
    return δ²β 
end

function _compute_δ³poly(t, A)
    δ³β = Expr(:call, :+, )
    for i in eachindex(A)
        δ³Aᵢ = (i-1) * (i-2) * (i-3) * A[i]
        if i > 3 && !(δ³Aᵢ ≈ 0.0)
            push!(
                δ³β.args, 
                (i-4) == 0 
                    ? δ³Aᵢ : Expr(:call, :*, δ³Aᵢ, Expr(:call, :^, t, Int(i-4)))
            )
        end    
    end
    if length(δ³β.args) == 1
        return :(0.0)
    end
    return δ³β 
end

function _compute_θ(t, Θm, n)
    Θᵢ = Expr(:call, :SVector, )
    lnz = []
    for i in 1:n
        θ = _compute_poly(t, Θm[i, 1:end])
        θ == 0.0 ? push!(lnz, 0) : push!(lnz, 1)
        push!(Θᵢ.args, θ)
    end
    return Θᵢ, lnz
end

function _compute_δθ(t, Θm, n)
    δΘᵢ = Expr(:call, :SVector, )
    lnz = []
    for i in 1:n
        θ = _compute_δpoly(t, Θm[i, 1:end])
        θ == 0.0 ? push!(lnz, 0) : push!(lnz, 1)
        push!(δΘᵢ.args, θ)
    end
    return δΘᵢ, lnz
end

function _compute_δ²θ(t, Θm, n)
    δ²Θᵢ = Expr(:call, :SVector, )
    lnz = []
    for i in 1:n
        θ = _compute_δ²poly(t, Θm[i, 1:end])
        θ == 0.0 ? push!(lnz, 0) : push!(lnz, 1)
        push!(δ²Θᵢ.args, θ)
    end
    return δ²Θᵢ, lnz
end

function _compute_δ³θ(t, Θm, n)
    δ³Θᵢ = Expr(:call, :SVector, )
    lnz = []
    for i in 1:n
        θ = _compute_δ³poly(t, Θm[i, 1:end])
        θ == 0.0 ? push!(lnz, 0) : push!(lnz, 1)
        push!(δ³Θᵢ.args, θ)
    end
    return δ³Θᵢ, lnz
end

function _add_sinusoidal!(β, B, f, i, vec)
    push!(
        β.args, 
        Expr(
            :call, 
            :*,
            B,
            Expr(:call, f, Expr(:call, :getindex, vec, i))
        )
    )
    nothing
end

function _iau_poly(tp, A, fp)
    βp = _compute_poly(tp, A)
    δβp = _compute_δpoly(tp, A)
    δ²βp = _compute_δ²poly(tp, A)
    δ³βp = _compute_δ³poly(tp, A)

    δβp_ = δβp == 0.0 ? :() : Expr(:call, :*, δβp, 1/fp)
    δ²βp_ = δ²βp == 0.0 ? :() : Expr(:call, :*, δ²βp, 1/fp^2)
    δ³βp_ = δ³βp == 0.0 ? :() : Expr(:call, :*, δ³βp, 1/fp^3)

    return βp, δβp_, δ²βp_, δ³βp_
end

function _iau_sine(ts, B, Θ, χ, fs)

    δχ = χ == :cos ? :sin : :cos
    δ²χ = χ == :cos ? :cos : :sin
    δ³χ = χ == :cos ? :sin : :cos

    sχ = χ == :cos ? -1.0 : 1.0
    s²χ = χ == :cos ? -1.0 : -1.0
    s³χ = χ == :cos ? 1.0 : -1.0

    βs = Expr(:call, :+, )
    δβs = Expr(:call, :+, )
    δ²βs = Expr(:call, :+, )
    δ³βs = Expr(:call, :+, )

    for i in eachindex(B)
        Bᵢ = B[i]
        if !(Bᵢ ≈ 0.0)
            _add_sinusoidal!(βs, Bᵢ, χ, i, :Θᵢ)

            # Θᵢ = Expr(:call, :getindex, :Θᵢ, i)
            # δΘᵢ = Expr(:call, :getindex, :δΘᵢ, i)
            # δ²Θᵢ = Expr(:call, :getindex, :δ²Θᵢ, i)
            # δ³Θᵢ = Expr(:call, :getindex, :δ³Θᵢ, i)

            # first derivative
            push!(
                δβs.args,
                Expr(
                    :call, 
                    :*, 
                    Bᵢ, 
                    Expr(:call, :getindex, :δΘᵢ, i), 
                    sχ, 
                    Expr(:call, δχ, Expr(:call, :getindex, :Θᵢ, i)) 
                )
            )

            # second derivative
            # δΘᵢ² = Expr(:call, :getindex, :δΘᵢ², i)
            push!(
                δ²βs.args, 
                Expr(
                    :call, 
                    :*,
                    Bᵢ,
                    Expr(
                        :call, 
                        :+, 
                        Expr(
                            :call, 
                            :*, 
                            sχ, 
                            Expr(:call, :getindex, :δ²Θᵢ, i), 
                            Expr(:call, δχ, Expr(:call, :getindex, :Θᵢ, i))),
                        Expr(
                            :call, 
                            :*, 
                            s²χ, 
                            Expr(:call, :getindex, :δΘᵢ², i), 
                            Expr(:call, δ²χ, Expr(:call, :getindex, :Θᵢ, i)))
                    )
                )
            )

            # third derivative
            # δΘᵢ³ = Expr(:call, :getindex, :δΘᵢ³, i)
            push!(
                δ³βs.args, 
                Expr(
                    :call, 
                    :*,
                    Bᵢ, 
                    Expr(
                        :call,
                        :+,
                        Expr(
                            :call, 
                            :*, 
                            sχ, 
                            Expr(:call, :getindex, :δ³Θᵢ, i), 
                            Expr(:call, δχ, 
                            Expr(:call, :getindex, :Θᵢ, i))),
                        Expr(
                            :call, 
                            :*, 
                            3.0, 
                            s²χ, 
                            Expr(:call, :getindex, :δΘᵢ, i), 
                            Expr(:call, :getindex, :δ²Θᵢ, i), 
                            Expr(:call, δ²χ, Expr(:call, :getindex, :Θᵢ, i))),
                        Expr(
                            :call, 
                            :*, 
                            s³χ, 
                            Expr(:call, :getindex, :δΘᵢ³, i), 
                            Expr(:call, δ³χ, 
                            Expr(:call, :getindex, :Θᵢ, i)))
                    )
                )
            )

        end
    end

    δβs_ = Expr(:call, :*, δβs, 1/fs)
    δ²βs_ = Expr(:call, :*, δ²βs, 1/fs^2)
    δ³βs_ = Expr(:call, :*, δ³βs, 1/fs^3)

    return βs, δβs_, δ²βs_, δ³βs_
end

function _iau_angles(tp::Symbol, ts::Symbol, A, B, Θ, χ, nuts, max_deg, fp=1.0, fs=1.0)
    # β = ∑ Aᵢ⋅tpⁱ + ∑ Bᵢ χ(∑ θᵢⱼ⋅tsʲ)

    β = Expr(:call, :+, )
    δβ = Expr(:call, :+, )
    δ²β = Expr(:call, :+, )
    δ³β = Expr(:call, :+, )
    
    n = length(Θ)÷(max_deg+1)
    Θm = reshape(Θ, (max_deg+1, n))'
    Θᵢ, _ = _compute_θ(ts, Θm, n)
    δΘᵢ, _ = _compute_δθ(ts, Θm, n)
    δ²Θᵢ, _ = _compute_δ²θ(ts, Θm, n)
    δ³Θᵢ, _ = _compute_δ³θ(ts, Θm, n) 

    # polynomial 
    if length(A) > 0

        βp, δβp, δ²βp, δ³βp = _iau_poly(tp, A, fp)

        if length(βp.args) > 1 push!(β.args, βp) end 
        if length(δβp.args) > 1 push!(δβ.args, δβp) end 
        if length(δ²βp.args) > 1 push!(δ²β.args, δ²βp) end 
        if length(δ³βp.args) > 1 push!(δ³β.args, δ³βp) end 

    end

    # sinusoidal
    if nuts && !(all(B .== 0.0)) 

        βs, δβs, δ²βs, δ³βs = _iau_sine(ts, B, Θ, χ, fs)

        if length(βs.args) > 1 push!(β.args, βs) end 
        if length(δβs.args) > 1 push!(δβ.args, δβs) end 
        if length(δ²βs.args) > 1 push!(δ²β.args, δ²βs) end 
        if length(δ³βs.args) > 1 push!(δ³β.args, δ³βs) end 

    end

    if length(β.args) ≤ 1
        β = :(0.0)
    end 
    if length(δβ.args) ≤ 1
        δβ = :(0.0)
    end 
    if length(δ²β.args) ≤ 1
        δ²β = :(0.0)
    end 
    if length(δ³β.args) ≤ 1
        δ³β = :(0.0)
    end 

    return β, δβ, δ²β, δ³β, Θᵢ, δΘᵢ, δ²Θᵢ, δ³Θᵢ

end

function _iau_angles(tp::Symbol, ts::Symbol,  p::PrecessionNutationComponent, fp=1.0, fs=1.0)
    return _iau_angles(tp, ts, p.A, p.B, p.Θ, p.fun, p.nuts, p.max_phase_deg, fp, fs)
end

function _iau_angles(p)
    return (
        _iau_angles(:T, :T, getproperty(p, :ra), Tempo.CENTURY2SEC, Tempo.CENTURY2SEC), 
        _iau_angles(:T, :T, getproperty(p, :dec), Tempo.CENTURY2SEC, Tempo.CENTURY2SEC),
        _iau_angles(:d, :T, getproperty(p, :pm), Tempo.DAY2SEC, Tempo.CENTURY2SEC)
    )
end

function orient_planets_angles(p, name)

    angles = _iau_angles(p)

    a, ad, add, addd, Θᵢ, δΘᵢ, δ²Θᵢ, δ³Θᵢ = angles[1]
    del, dd, ddd, dddd, _, _, _, _ = angles[2]
    w, wd, wdd, wddd, _, _, _, _ = angles[3]

    bname = lowercase(String(name))

    fname = Symbol("orient_angles_$bname")
    fdname = Symbol("orient_d_angles_$bname")
    fddname = Symbol("orient_dd_angles_$bname")

    return eval(quote 
        function ($fname)(sec::Number)
            d = sec/Tempo.DAY2SEC
            T = d/Tempo.CENTURY2DAY
            Θᵢ = $Θᵢ
            return $a, $del, $w
        end,
        function ($fdname)(sec::Number)
            d = sec/Tempo.DAY2SEC
            T = d/Tempo.CENTURY2DAY
            Θᵢ = $(Θᵢ)
            δΘᵢ = $(δΘᵢ)
            return $a, $del, $w, $ad, $dd, $wd
        end,
        function ($fddname)(sec::Number)
            d = sec/Tempo.DAY2SEC
            T = d/Tempo.CENTURY2DAY
            Θᵢ = $(Θᵢ)
            δΘᵢ = $(δΘᵢ)
            δ²Θᵢ = $(δ²Θᵢ)
            δΘᵢ² = δΘᵢ .* δΘᵢ
            return $a, $del, $w, $ad, $dd, $wd, $add, $ddd, $wdd
        end
    end)
end