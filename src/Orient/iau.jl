function parse_anglestr(A::NV, B::Union{Nothing, NV}, Θ::Union{Nothing, NM}, 
    t::Symbol, χ::Symbol; conv::Real=π/180) where {NV<:AbstractArray, NM<:AbstractArray}
    # β = ∑AᵢTⁱ + ∑ Bᵢ χ(θ₀ᵢ + θ₁ᵢ t)

    if (B !== nothing && Θ !== nothing)
        if length(B) != size(Θ)[1]
            throw(error("[Orient] Rotation/nutation angles and coeffs" 
                * "must be of the same length."))
        end
        nuts = true
    end

    δχ = χ == :cos ? :sin : :cos
    sχ = χ == :cos ? -1.0 : 1.0

    # angle
    β = ""
    for (i, ai) in enumerate(A)
        Aᵢ = ai*conv
        if !(Aᵢ ≈ 0.0)
            β *= "+$Aᵢ"*"*$t"^(i-1)
        end
    end

    # rate 
    δβ = ""
    for (i, ai) in enumerate(A[2:end])
        δAᵢ = i*ai*conv
        if !(δAᵢ ≈ 0.0)
            δβ *= "+$(δAᵢ)"*"*$t"^(i-1)
        end
    end

    # nutation-precession angles 
    if nuts
        for (i, bi) in enumerate(B)
            Bᵢ = bi*conv
            if !(Bᵢ ≈ 0.0)
                θᵢ = "$(Θ[i, 1])+$(Θ[i, 2])*$t"
                β *= "+$(Bᵢ)*$χ($θᵢ)"
                δBᵢ = Bᵢ*sχ*Θ[i, 2]
                if !(δBᵢ ≈ 0.0)
                    δβ *= "+$(δBᵢ)*$(δχ)($θᵢ)"
                end
            end
            
        end
    end

    strip(β) == "" ? throw(error("[Orient] Angle cannot be parsed!" * 
        "There is nothing to parse!")) : ()
    strip(δβ) == "" ? δβ = "0.0" : ()
    
    β, δβ
end