"""
    fw_angles(::IAU2006, j2000ttc::N)

Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).

### Inputs
- IAU Model type
- `j2000ttc`  -- `TT` centuries since J2000

### Outputs
- `γ` -- F-W 1st angle -- `rad`
- `ϕ` -- F-W 2nd angle -- `rad`
- `ψ` -- F-W 3rd angle -- `rad`   
- `ε` -- F-W 4th angle -- `rad`

### References 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pfw06.c) library
"""
function fw_angles(::IAU2006, j2000ttc::N) where {N<:Number}
    γ = @evalpoly(
        j2000ttc,
        -0.052928,
        10.556378,
        0.4932044,
        -0.00031238,
        -0.000002788,
        0.0000000260,
    ) |> arcsec2rad

    ϕ = @evalpoly(
        j2000ttc,
        84381.412819,
        -46.811016,
        0.0511268,
        0.00053289,
        -0.000000440,
        -0.0000000176,
    ) |> arcsec2rad

    ψ = @evalpoly(
        j2000ttc,
        -0.041775,
        5038.481484,
        1.5584175,
        -0.00018522,
        -0.000026452,
        -0.0000000148,
    ) |> arcsec2rad

    ε = mean_obliquity(iau2006, j2000ttc)
    return γ, ϕ, ψ, ε

end

"""
    fw_matrix(γ, ϕ, ψ, ε)

Form rotation matrix given the Fukushima-Williams angles.

### References
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pmat06.c) software library
"""
function fw_matrix(γ, ϕ, ψ, ε)
    return angle_to_dcm(γ, :Z) * angle_to_dcm(ϕ, -ψ, -ε, :XZX)
end

"""
    fw_precession_bias(::IAU2006, j2000ttc::N) where {N<:Number}

Precession matrix (including frame bias) from GCRS to a specified
date, IAU 2006 model. This uses the Fukushima-Williams model.

### Inputs 
- `j2000ttc`-- `TT` centuries since J2000

### Output 
Rotation matrix. 

The matrix operates in the sense V(date) = RBP * V(GCRS), where the vector 
V(GCRS) is with respect to the Geocentric Celestial Reference System (IAU, 2000) 
and the vector V(date) is with respect to the mean equatorial triad of the 
given date.

### References:

- Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

- IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
    24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6. (2000)

- Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981

- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pmat06.c) software library
"""
function fw_precession_bias(q::M, j2000ttc::N) where {N<:Number, M<:IAU2006Model}
    # Bias-precession Fukushima-Williams angles.
    γ, ϕ, ψ, ε = fw_angles(q, j2000ttc)
    # form the matrix
    return fw_matrix(γ, ϕ, ψ, ε)
end

const ICRF2J2000_BIAS= fw_precession_bias(iau2006a, 0.0)

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