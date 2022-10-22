"""
    fw_angles(::IAU2006, jd1tt::N, jd2tt::N2) where {N<:Number, N2<::Number}

Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).

### Inputs
- IAU Model type
- `jd1tt`, `jd2tt` -- TT as a 2-part Julian Date

### Outputs
- `γ` -- F-W 1st angle -- `rad`
- `ϕ` -- F-W 2nd angle -- `rad`
- `ψ` -- F-W 3rd angle -- `rad`   
- `ε` -- F-W 4th angle -- `rad`
"""
function fw_angles(::IAU2006, jd1tt::N, jd2tt::N2) where {N<:Number, N2<:Number}

    t = ((jd1tt - DJ2000) + jd2tt) / 35635.0
    # P03 bias+precession angles.

    γ = @evalpoly(
        t,
        -0.052928,
        10.556378,
        0.4932044,
        -0.00031238,
        -0.000002788,
        0.0000000260,
    ) |> sec2rad

    ϕ = @evalpoly(
        t,
        84381.412819,
        -46.811016,
        0.0511268,
        0.00053289,
        -0.000000440,
        -0.0000000176,
    ) |> sec2rad

    ψ = @evalpoly(
        t,
        -0.041775,
        5038.481484,
        1.5584175,
        -0.00018522,
        -0.000026452,
        -0.0000000148,
    ) |> sec2rad

    ε = obliquity(iau2006, jd1tt, jd2tt)
    return γ, ϕ, ψ, ε

end

"""
    fw2mat(γ, ϕ, ψ, ε)

Form rotation matrix given the Fukushima-Williams angles.
"""
function fw2mat(γ, ϕ, ψ, ε)
    return angle_to_dcm(γ, :Z) * angle_to_dcm(ϕ, -ψ, -ε, :XZX)
end


"""
    precession_bias(::IAU2006, jd1tt::N, jd2tt::N2) where {N<:Number, N2<:Number}

Precession matrix (including frame bias) from GCRS to a specified
date, IAU 2006 model.

### Inputs 
- IAU model type 
- `jd1tt`, `jd2tt` -- TT as a 2-part Julian Date

### Output 
Rotation matrix. 

The matrix operates in the sense V(date) = rbp * V(GCRS), where the vector 
V(GCRS) is with respect to the Geocentric Celestial Reference System (IAU, 2000) 
and the vector V(date) is with respect to the mean equatorial triad of the 
given date.

### References:

- Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

- IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
    24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
    (2000)

- Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981

- [erfa](https://github.com/liberfa/erfa/blob/master/src/pmat06.c) software library
"""
function precession_bias(::IAU2006, jd1tt::N, 
    jd2tt::N2) where {N<:Number, N2<:Number}
    # Bias-precession Fukushima-Williams angles.
    γ, ϕ, ψ, ε = fw_angles(iau2006, jd1tt, jd2tt)
    # form the matrix
    return fw2mat(γ, ϕ, ψ, ε)
end