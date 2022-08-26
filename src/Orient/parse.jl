export orient_iau

function parse_angle_string(a, b, nuts::Bool, Θ, trig::Symbol, trigrate::Symbol, sign)
    # parse function
    D2R = π/180
    c = ""
    for (i, ri) in enumerate(a) c *= "+$(D2R*ri)"*"*T"^(i-1) end
    nuts && for i in range(1, length(b)) c*= "+$(D2R*b[i])*$trig($(D2R*Θ[i, 1])+$(D2R*Θ[i, 2])*T)" end 

    # parse function time derivative
    cd = ""
    for (i, ri) in enumerate(a[2:end]) 
        if !(ri*i ≈ 0.0)
            cd *= "+$(D2R*ri*i)"*"*T"^(i-1) 
        end 
    end
    nuts && for i in range(1, length(b)) cd*= "+$(D2R*Θ[i, 2]*sign*b[i])*$trigrate($(D2R*Θ[i, 1])+$(D2R*Θ[i, 2])*T)" end 
    if cd == "" cd = "0.0" end 
    return c, cd
end

"""
    orient_iau(pck::Dict, id2name::Dict)
    
This function automatically parse the body right ascention, declination and rotation 
and their rates from a parsed PCK2 file.

!!!! note 
    PCK *type 2* kernels are a `Basic.jl` modern rewriting of classical NAIF PCK kernels.
    Be sure to load the kernels of this format to parse the data.

!!!! note 
    This function requires the bodies associated types to be already parsed.
"""
function orient_iau(pck::Dict, id2name::Dict)
    @info "[Orient] Populating IAU constants for celestial bodies body-fixed frames."
    sv = Vector{Expr}()
    for id in keys(id2name)
        body = id2name[id]

        # RIGHT ASCENSION
        @debug "[Orient] Parsing $id body right-ascention."
        a = parse_pck_property("POLE_RA", id, pck)
        nuts = false
        Θ = nothing
        b = parse_pck_property("NUT_PREC_RA", id, pck)
        if !isnothing(b)
            Θ = parse_pck_property("NUT_PREC_RA_ANGLES", id, pck)
            isnothing(Θ) && throw(error("Impossible to parse `NUT_PREC_RA_ANGLES` for $id. Please load TPC kernel of type 2."))
            nuts = true
            Θ = transpose(reshape(Θ, (2, length(b))))
            @debug "[Orient] Body $id right-ascention has nutation and precession."
        end
        c, cd = parse_angle_string(a, b, nuts, Θ, :sin, :cos, 1.0)
        
        f = quote
            function body_right_ascention(::$body, T::Float64) 
                $(Meta.parse(c))
            end
            function body_right_ascention_rate(::$body, T::Float64) 
                $(Meta.parse(cd))
            end
        end
        push!(sv, f)

        # DECLINATION
        @debug "[Orient] Parsing $id body declination."
        a = parse_pck_property("POLE_DEC", id, pck)
        nuts = false
        Θ = nothing
        b = parse_pck_property("NUT_PREC_DEC", id, pck)
        if !isnothing(b)
            Θ = parse_pck_property("NUT_PREC_DEC_ANGLES", id, pck)
            isnothing(Θ) && throw(error("Impossible to parse `NUT_PREC_DEC_ANGLES` for $id. Please load TPC kernel of type 2."))
            nuts = true
            Θ = transpose(reshape(Θ, (2, length(b))))
            @debug "[Orient] Body $id declination has nutation and precession."
        end
        c, cd = parse_angle_string(a, b, nuts, Θ, :cos, :sin, -1.0)
        
        f = quote 
            function body_declination(::$body, T::Float64) 
                $(Meta.parse(c))
            end
            function body_declination_rate(::$body, T::Float64) 
                $(Meta.parse(cd))
            end
        end
        push!(sv, f)

        # ROTATION
        @debug "[Orient] Parsing $id body rotation."
        a = parse_pck_property("PM", id, pck)
        nuts = false
        Θ = nothing
        b = parse_pck_property("NUT_PREC_PM", id, pck)
        if !isnothing(b)
            Θ = parse_pck_property("NUT_PREC_PM_ANGLES", id, pck)
            isnothing(Θ) && throw(error("Impossible to parse `NUT_PREC_PM_ANGLES` for $id. Please load TPC kernel of type 2."))
            nuts = true
            Θ = transpose(reshape(Θ, (2, length(b))))
            @debug "[Orient] Body $id rotation has nutation and precession."
        end

        c, cd = parse_angle_string(a, b, nuts, Θ, :sin, :cos, 1.0)
        f = quote
            function body_rotation_angle(::$body, T::Float64) 
                $(Meta.parse(c))
            end
            function body_rotation_rate(::$body, T::Float64) 
                $(Meta.parse(cd))
            end
        end
        push!(sv, f)

    end
    @info "[Orient] IAU α, δ, ω populated for $(values(id2name))."
    sv
end
