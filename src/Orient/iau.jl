export generate_iauangles!, parse_iauconstants

import Basic.Utils: genf_psnginfst

function parse_iauanglestr(A::NV, B::Union{Nothing, NV}, Θ::Union{Nothing, NM}, 
    t::Symbol, χ::Symbol; conv::AbstractFloat=π/180) where {NV<:AbstractArray, NM<:AbstractArray}
    # β = ∑AᵢTⁱ + ∑ Bᵢ χ(θ₀ᵢ + θ₁ᵢ t)
    nuts = false
    if (B !== nothing && Θ !== nothing)
        if length(B) > size(Θ)[1]
            # the length of B and Θ are NOT equal because, according to IAU
            throw(error("[Orient] Rotation/nutation angles and coeffs" 
                * "must compatible in size as for IAU standards."))
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

function parse_iauconstants(file::String)
    tpc = load(TPC(file))
    parse_iauconstants(collect(keys(tpc)), tpc)
end

function parse_iauconstants(bodiesid::Vector{N}, 
    data::D) where {N <: Integer, D <: AbstractDict}

    parsed = Dict{N,Union{Dict{String,Dict{Symbol,Union{Array,Nothing,Symbol}}},Nothing}}()
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

        !(isempty(di)) ? push!(parsed, bid => di) : push!(parsed, bid => nothing)
    end
    parsed
end

function generate_iauangles!(gen::String, bid::N, 
    iaudata::D; conv::N2=π/180) where {N<:Integer, D<:AbstractDict, N2<:AbstractFloat}
    @debug "[Orient] Generating IAU model for $bid"

    if haskey(iaudata, bid)
        gen *= "\n"
        gen *= "#%ORIENT::$bid\n"
        iauconst = iaudata[bid]
        if !(iauconst === nothing)
            for (angle, data) in pairs(iauconst)
                gen *= "#%ORIENT::$bid/$angle\n"
                β, δβ = parse_iauanglestr(data[:A], data[:B], data[:Θ], :T, data[:χ]; conv=conv)
                angle = "orient_$angle"
                gen *= genf_psnginfst(:Orient, angle == "orient_rotation" ? "orient_rotation_angle" : angle, β, 
                    (nothing, Val{bid}), (:T, :N); wherestr="N<:AbstractFloat")
                gen *= genf_psnginfst(:Orient, join((angle,"rate"),"_"), δβ, 
                    (nothing, Val{bid}), (:T, :N); wherestr="N<:AbstractFloat")
            end 
        else 
            for fun in ("orient_declination", "orient_declination_rate", 
                "orient_right_ascension", "orient_right_ascension_rate", 
                "orient_rotation_angle", "orient_rotation_rate")
                errorprop = join(uppercasefirst.(split(fun,"_")[2:end]), " ")
                gen *= genf_psnginfst(:Orient, fun, 
                    """throw(error("[Orient] IAU `$(errorprop)` cannot be computed for $bid."))""", 
                    (nothing, Val{bid}), (:T, :N); wherestr="N<:AbstractFloat")

            end
        end

    end
    gen
end 
