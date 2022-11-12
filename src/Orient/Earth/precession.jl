export orient_precession_bias

"""
    fw_angles(::IAU2006Model, t::N) where {N<:Number}

Compute the precession angles in radians, following the IAU 2006 
Fukushima-Williams 4-angle formulation.

### Inputs
- IAU Model type
- `t`  -- Terrestrial Time `TT` Julian centuries since J2000

### Outputs
- `γ` -- F-W 1st angle
- `ϕ` -- F-W 2nd angle
- `ψ` -- F-W 3rd angle   
- `ε` -- F-W 4th angle

### References 
- Luzum, B. and Petit G. (2012), _The IERS Conventions (2010)_, 
[IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Wallace P. T. and Capitaine N. (2006), _Precession-nutation procedures consistent with 
IAU 2006 resolutions_, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pfw06.c) library
"""
function fw_angles(::IAU2006Model, t::N) where {N<:Number}
    
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

Form the Nutation-Precession-Bias (NPB) rotation matrix given the 
Fukushima-Williams angles, expressed in radians.

### Notes 
The present function can construct three different matrices depending on which 
angles are supplied as arguments: 

- To obtain the NPB matrix, generate the four standard FW precession angles 
(̄γ, ̄ϕ, ̄ψ, ϵₐ) then generate the nutation components Δψ and Δϵ and add them to 
̄ψ, ϵₐ. Finally, call the present functions using those four angles as arguments. 

- To obtain the precession-frame bias matrix (PB), generate the four standard FW 
precession angles and call the present function. 

- To obtain the frame bias matrix (B), generate the four standard FW precession angles 
at date J2000.0 and call this function.

The remaining nutation-only and precession only matrices can be obtained by 
combining these three appropriately. 

### References
- Wallace P. T. and Capitaine N. (2006), _Precession-nutation procedures consistent with 
IAU 2006 resolutions_, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pmat06.c) library
"""
function fw_matrix(γ, ϕ, ψ, ϵ)
    angle_to_dcm(-ϵ, :X)*angle_to_dcm(γ, ϕ, -ψ, :ZXZ)
end


"""
    orient_precession_bias(q::IAU2006, t::N) where {N<:Number}

Form the precession-frame bias (PB) matrix from GCRS to a specified
date, using the Fukushima-Williams angles and following the IAU 2006 model.

### Inputs 
- `t`-- Terrestrial Time `TT` Julian centuries since J2000

### Notes
The matrix operates in the sense V(date) = RBP * V(GCRS), where the vector 
V(GCRS) is with respect to the Geocentric Celestial Reference System (IAU, 2000) 
and the vector V(date) is with respect to the mean equatorial triad of the 
given date.

### References:
- IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
    24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6. (2000)
- Wallace P. T. and Capitaine N. (2006), _Precession-nutation procedures consistent with 
    IAU 2006 resolutions_, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pmat06.c) software library
"""
function orient_precession_bias(q::M, t::N) where {N<:Number, M<:IAU2006Model}
    # Bias-precession Fukushima-Williams angles.
    γ, ϕ, ψ, ε = fw_angles(q, t)
    # form the matrix
    return fw_matrix(γ, ϕ, ψ, ε)
end

const ICRF2J2000_BIAS = orient_precession_bias(iau2006a, 0.0)


"""
    ICRF2J2000_BIAS

Rotation matrix for the rotation from the International Celestial Reference Frame 
(`ICRF`) and the Mean Dynamical Equator and Equinox at J2000.0 (`MEME2000`).

### References
- Hilton, James L., and Catherine Y. Hohenkerk. -- Rotation matrix from the mean 
    dynamical equator and equinox at J2000. 0 to the ICRS. -- Astronomy & Astrophysics 
    413.2 (2004): 765-770. DOI: [10.1051/0004-6361:20031552](https://www.aanda.org/articles/aa/pdf/2004/02/aa3851.pdf)
- [SOFA docs](https://www.iausofa.org/2021_0512_C/sofa/sofa_pn_c.pdf)
"""
ICRF2J2000_BIAS