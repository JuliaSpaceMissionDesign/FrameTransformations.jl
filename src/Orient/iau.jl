function parse_iauanglestr(A::NV, B::Union{Nothing, NV}, Θ::Union{Nothing, NM}, 
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

    β = ""
    δβ = ""

    if A !== nothing
        for (i, ai) in enumerate(A)
            Aᵢ = ai*conv
            if !(Aᵢ ≈ 0.0)
                β *= "+$Aᵢ"*"*$t"^(i-1)
            end
        end
        for (i, ai) in enumerate(A[2:end])
            δAᵢ = i*ai*conv
            if !(δAᵢ ≈ 0.0)
                δβ *= "+$(δAᵢ)"*"*$t"^(i-1)
            end
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

function _format_iauconst(hasnuts, dbid, nuts, prop, χ)
    pole = prop != :pm ? Symbol("pole_$prop") : :pm
    nutprec = Symbol("nut_prec_$prop")
    Dict{Symbol, Union{Array, Nothing, Symbol}}(
        :A => haskey(dbid, pole) ? dbid[pole] : nothing,
        :B => hasnuts && haskey(dbid, nutprec) ? dbid[nutprec] : nothing,
        :Θ => nuts,
        :χ => χ
    )
end

function parse_iauconstants(bodiesid::Vector{N}, 
    data::D) where {N <: Integer, D <: AbstractDict}

    parsed = Dict{N, Dict{String, Dict{Symbol, Union{Array, Nothing, Symbol}}}}()
    for bid in bodiesid
        @debug "[Orient] Parsing constants for object with NAIF ID: $bid."
        sid = "$(bid)"
        nutsid = bid
        if bid < 1000 && bid > 100
            nutsid = parse(N, sid[1])
        end
        
        di = Dict{String, Dict{Symbol, Union{Array, Nothing, Symbol}}}()
        if bid < 9
           @warn "[Orient] Cannot orient IAU frames for $bid - IGNORED."
        else
            dbid = data[bid]

            # nutation precession angles 
            nuts = haskey(data[nutsid], :nut_prec_angles) ? data[nutsid][:nut_prec_angles] : nothing
            hasnuts = !(nuts === nothing)
            nuts = hasnuts ? collect(transpose(reshape(nuts, (2, Int(length(nuts)/2))))) : nothing
            
            # right-ascension
            push!(di, 
                "right_ascension" => _format_iauconst(hasnuts, 
                    dbid, nuts, :ra, :sin))

            # declination
            push!(di, 
                "declination" => _format_iauconst(hasnuts, 
                    dbid, nuts, :dec, :cos))

            # polar motion
            push!(di, 
                "rotation" => _format_iauconst(hasnuts, 
                    dbid, nuts, :pm, :sin))
            # adjust rotation from day to centuries
            if !(di["rotation"][:A] === nothing) && length(di["rotation"][:A]) > 1
                @. di["rotation"][:A][2:end] *= 36525.0
            end 
        end

        !(isempty(di)) ? push!(parsed, bid => di) : ()
    end
    parsed
end