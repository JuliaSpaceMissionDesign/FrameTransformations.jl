export orient_bias_precession, 
       orient_bpn


"""
    fw_angles(m::IAU2006Model, t::Number) 

Compute the precession angles in radians, following the IAU 2006 Fukushima-Williams 4-angle 
formulation at time `t` expressed in `TT` Julian centuries since [`J2000`](@ref).

### Outputs
- `γ` -- F-W 1st angle
- `ϕ` -- F-W 2nd angle
- `ψ` -- F-W 3rd angle   
- `ε` -- F-W 4th angle

### References 
- Luzum, B. and Petit G. (2012), The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pfw06.c) library
"""
function fw_angles(::IAU2006Model, t::Number) 
    
    γ = @evalpoly(
        t,
        -0.052928,
        10.556378,
         0.4932044,
        -0.00031238,
        -0.000002788,
         0.0000000260,
    ) |> arcsec2rad

    ϕ = @evalpoly(
        t,
        84381.412819,
          -46.811016,
            0.0511268,
            0.00053289,
           -0.000000440,
           -0.0000000176,
    ) |> arcsec2rad

    ψ = @evalpoly(
        t,
          -0.041775,
        5038.481484,
           1.5584175,
          -0.00018522,
          -0.000026452,
          -0.0000000148,
    ) |> arcsec2rad

    ϵ = orient_obliquity(iau2006a, t)
    return γ, ϕ, ψ, ϵ

end


"""
    fw_matrix(γ, ϕ, ψ, ε)

Form the Nutation-Precession-Bias (NPB) rotation matrix given the Fukushima-Williams angles, 
expressed in radians.

The present function can construct three different matrices depending on which angles are 
supplied as arguments: 

- **NPB**: To obtain the Nutation-Precession-Bias matrix, generate the four standard FW 
    precession angles (̄γ, ̄ϕ, ̄ψ, ϵₐ) then generate the nutation components Δψ and Δϵ and add them 
    to ̄ψ, ϵₐ. Finally, call the present functions using those four angles as arguments. 

- **PB**: To obtain the precession-frame bias matrix, generate the four standard FW precession 
    angles and call the present function. 

- **B**: To obtain the frame bias matrix, generate the four standard FW precession angles at 
    date J2000.0 and call this function.

The remaining nutation-only and precession only matrices can be obtained by combining these 
three appropriately. 

### References
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/fw2m.c) library
"""
function fw_matrix(γ, ϕ, ψ, ϵ)
    angle_to_dcm(-ϵ, :X)*angle_to_dcm(γ, ϕ, -ψ, :ZXZ)
end


"""
    precession_angles(m::IAU1980Model, t::Number)

Compute the precession angles from Lieske et al. 1977 model, in radians, at time `t` 
expressed in `TT` Julian centuries since [`J2000`](@ref).

### References 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/bp00.c) software library
"""
function precession_angles(::IAU1980Model, t::Number)
    # Compute Precession Angles from Lieske et al. 1977 
    ψₐ = @evalpoly(t, 0.0, 5038.7784, -1.07259, -0.001147) |> arcsec2rad
    ωₐ = @evalpoly(t, 84381.448, 0.0, 0.05127, -0.007726) |> arcsec2rad
    χₐ = @evalpoly(t, 0.0, 10.5526, -2.38064, -0.001125) |> arcsec2rad

    return ψₐ, ωₐ, χₐ
end


"""
    precession_rate(m::IAU2000Model, t::Number)

Compute the precession-rate part of the IAU 2000 precession-nutation models, in radians, at 
time `t` expressed as `TT` Julian centuries since [`J2000`](@ref).

### References 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pr00.c) software library
"""
function precession_rate(::IAU2000Model, t::Number)

    Δψₚ = arcsec2rad(-0.29965t)
    Δϵₚ = arcsec2rad(-0.02524t)

    return Δψₚ, Δϵₚ
end


"""
    frame_bias(::IAU2000Model)

Compute the frame bias components of the IAU 2000 precession-nutation models, in radians.

### Notes 
- The frame bias corrections in longitude and obliquity are required to correct for the 
  offset between the GCRS pole and the mean [`J2000`](@ref) pole. They define, with respect 
  to the GCRS axes, a J2000 mean pole that is consistent with teh IAU 2000A precession-nutation 
  model. 

- The function also returns an offset in right ascension taken from Chapront et al. (2002),
  necessary to completely describe the frame bias, but that is not part of the original IAU 
  2000A model.

### References 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/bi00.c) software library
"""
function frame_bias(::IAU2000Model)

    # The frame bias corrections in longitude and obliquity
    Δψb = arcsec2rad(-0.041775) 
    Δϵb = arcsec2rad(-0.0068192)

    # The ICRS RA of the J2000.0 equinox (Chapront et al., 2002) 
    Δα₀ = arcsec2rad(-0.0146)
    return Δψb, Δϵb, Δα₀
end


"""
    orient_bias_precession(m::IAUModel, t::Number)

Form the precession-frame bias (PB) matrix that transforms vectors from the GCRS to the 
mean of date, following the IAU Model `m` at time `t` expressed as `TT` Julian centuries 
since [`J2000`](@ref).

The function has been implemented for the `IAU2000` and `IAU2006` models.

### References:
- IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc. 24th General Assembly, 
  Manchester, UK.  Resolutions B1.3, B1.6. (2000)
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with IAU 
  2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- ERFA [pmat06](https://github.com/liberfa/erfa/blob/master/src/pmat06.c) function.
- ERFA [pmat00](https://github.com/liberfa/erfa/blob/master/src/pmat00.c) function.
"""
function orient_bias_precession(m::IAU2006Model, t::Number)
    # Bias-precession Fukushima-Williams angles.
    γ, ϕ, ψ, ε = fw_angles(m, t)

    # Form the matrix
    return fw_matrix(γ, ϕ, ψ, ε)
end

function orient_bias_precession(m::IAU2000Model, t::Number)

    # J2000.0 obliquity (Lieske et al. 1977)
    ϵ₀ = arcsec2rad(84381.448); # arcseconds 

    # Frame bias matrix: it transforms vectors from GCRS to mean J2000.0
    δψᵦ, δϵᵦ, δα₀ = frame_bias(m) 
    Rᵦ = angle_to_dcm(δα₀, δψᵦ*sin(ϵ₀), -δϵᵦ, :ZYX)

    # Precession angles
    ψₐ, ωₐ, χₐ = precession_angles(iau1980, t)

    # Apply IAU 2000 precession corrections
    Δψₚ, Δϵₚ = precession_rate(m, t)

    ψₐ += Δψₚ
    ωₐ += Δϵₚ

    # Precession matrix: it transforms from mean J2000.0 to mean of date
    Rₚ = angle_to_dcm(χₐ, :Z)*angle_to_dcm(ϵ₀, -ψₐ, -ωₐ, :XZX)

    return Rₚ*Rᵦ
end


""" 
    orient_bias_precession_nutation(m::IAUModel, t::Number)

Compute the equinox-based bias-precession-nutation matrix using the IAU Model `m` procedures 
at time `t` expressed in `TT` Julian centuries since [`J2000`](@ref).

The function has been implemented for the `IAU2000` and `IAU2006` models.

!!! note 
    The computed matrix rotates a vector from the GCRS to the true equatorial triad of date.

### References 
- ERFA [pn06](https://github.com/liberfa/erfa/blob/master/src/pn06.c) function
- ERFA [pn00](https://github.com/liberfa/erfa/blob/master/src/pn00.c) function
"""
function orient_bias_precession_nutation(m::IAU2006Model, t::Number)
    # Computes Fukushima-Williams angles
    γ, ϕ, ψ, ϵ = fw_angles(m, t)

    # Computes IAU 2000 nutation components 
    Δψ, Δϵ = orient_nutation(m, t) 

    # Equinox-based Bias-precession-nutation matrix, with 
    # applied IAU-2006 compatible nutations 
    fw_matrix(γ, ϕ, ψ + Δψ, ϵ + Δϵ)    
end

function orient_bias_precession_nutation(m::IAU2000Model, t::Number)
    _, Δϵₚ = precession_rate(m, t)
    ϵₐ = orient_obliquity(iau1980, t) + Δϵₚ

    # Frame bias and precession 
    RPB = orient_bias_precession(m, t)

    # Nutation Matrix
    Δψ, Δϵ = orient_nutation(m, t)
    RN = angle_to_dcm(ϵₐ, -Δψ, -(ϵₐ+Δϵ), :XZX)

    # Bias-precession-nutation matrix! 
    return RN*RPB    
end

