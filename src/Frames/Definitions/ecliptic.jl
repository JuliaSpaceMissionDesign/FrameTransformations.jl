export add_axes_eclipj2000!, add_axes_meme2000!

# --------------------------------------------------------
# DCMs
# --------------------------------------------------------

"""
    DCM_ICRF2J2000_BIAS

DCM for the rotation from the International Celestial Reference Frame (`ICRF`) and the 
Mean Dynamical Equator and Equinox of J2000.0 (`MEME2000`).

### References
- Hilton, James L., and Catherine Y. Hohenkerk. -- Rotation matrix from the mean 
    dynamical equator and equinox at J2000. 0 to the ICRS. -- Astronomy & Astrophysics 
    413.2 (2004): 765-770. DOI: [10.1051/0004-6361:20031552](https://www.aanda.org/articles/aa/pdf/2004/02/aa3851.pdf)
- [SOFA docs](https://www.iausofa.org/2021_0512_C/sofa/sofa_pn_c.pdf)
"""
const DCM_ICRF2J2000_BIAS = Orient.orient_precession_bias(Orient.iau2006a, 0.0)

"""
    DCM_J20002ECLIPJ2000

DCM for the rotation from the Mean Dynamical Equator of J2000 (`MEME2000`) to the 
Mean Ecliptic Equinox. This corresponds to the transformation `J2000 -> ECLIPJ2000` 
in the SPICE toolkit.
"""
const DCM_J20002ECLIPJ2000 = angle_to_dcm(
    Orient.orient_obliquity(Orient.iau2006a, 0.0), :X
)

"""
    DCM_ICRF2ECLIPJ2000

DCM for the rotation from the International Celestial Reference Frame (`ICRF`) to the 
Mean Ecliptic Equinox of J2000 (`ECLIPJ2000`).
"""
const DCM_ICRF2ECLIPJ2000 = DCM_ICRF2J2000_BIAS * DCM_J20002ECLIPJ2000

# --------------------------------------------------------
# TRANSFORMATIONS
# --------------------------------------------------------

function orient_axes_icrf_to_meme2000(t::Number)
    return DCM_ICRF2J2000_BIAS
end

function orient_axes_d_icrf_to_meme2000(t::Number)
    return DCM_ICRF2J2000_BIAS, DCM(0.0I)
end

function orient_axes_dd_icrf_to_meme2000(t::Number)
    return DCM_ICRF2J2000_BIAS, DCM(0.0I), DCM(0.0I)
end

function orient_axes_icrf_to_eclipj2000(t::Number)
    return DCM_ICRF2ECLIPJ2000
end

function orient_axes_d_icrf_to_eclipj2000(t::Number)
    return DCM_ICRF2ECLIPJ2000, DCM(0.0I)
end

function orient_axes_dd_icrf_to_eclipj2000(t::Number)
    return DCM_ICRF2ECLIPJ2000, DCM(0.0I), DCM(0.0I)
end

function orient_axes_icrf_to_mememod(sec::Number; model::Orient.IAU2006Model=iau2006b)
    γ, ϕ, ψ, ε = Orient.fw_angles(model, sec)
    R = Orient.fw_matrix(γ, ϕ, ψ, ε)
    return Rotation(R)
end

function orient_axes_meme2000_to_eclipj2000(t::Number)
    return DCM_J20002ECLIPJ2000
end

function orient_axes_d_meme2000_to_eclipj2000(t::Number)
    return DCM_J20002ECLIPJ2000, DCM(0.0I)
end

function orient_axes_dd_meme2000_to_eclipj2000(t::Number)
    return DCM_J20002ECLIPJ2000, DCM(0.0I), DCM(0.0I)
end

function add_axes_meme2000!(frames, axes, parent)
    pname = axes_name(parent)

    if pname != :ICRF 
        throw(
            ErrorException("Mean Equator, Mean Equinox of J2000 (MEME2000) axes could be defined only" * 
            " w.r.t the International Celestial Reference Frame(ICRF)")
        )
    end

    add_axes_inertial!(frames, axes; parent=parent, dcm=DCM_ICRF2J2000_BIAS)
end

function add_axes_eclipj2000!(frames, axes, parent)
    pname = axes_name(parent)

    if pname == :ICRF 
        dcm = DCM_J20002ECLIPJ2000 * DCM_ICRF2J2000_BIAS 
    elseif pname == :MEME2000 
        dcm = DCM_J20002ECLIPJ2000
    else
        throw(
            ErrorException("Ecliptic Equinox of J2000 (ECLIPJ2000) axes could not be defined" *
            " w.r.t $pname axes. `ICRF` or `MEME2000` could be instead used.")
        )
    end
    add_axes_inertial!(frames, axes; parent=parent, dcm=dcm)
end
