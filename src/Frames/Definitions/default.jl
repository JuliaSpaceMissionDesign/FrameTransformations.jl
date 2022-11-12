# ------------------------------------------------------------------------------
#                               ICRF <-> MEME2000
# ------------------------------------------------------------------------------

function Rotation(
    origin::InternationalCelestialReferenceFrame, 
    target::MeanEquatorMeanEquinoxJ2000, e::Epoch)
    return Rotation(origin, target, Orient.ROTMAT_ICRF2J2000_BIAS)
end

function Rotation(
    origin::MeanEquatorMeanEquinoxJ2000, 
    target::InternationalCelestialReferenceFrame, e::Epoch)
    Rotation(origin, target, Orient.ROTMAT_ICRF2J2000_BIAS')
end

connect!(ICRF, MEME2000)

const ROTMAT_ICRF2J2000_BIAS = Orient.orient_precession_bias(Orient.iau2006a, 0.0)

"""
    ROTMAT_ICRF2J2000_BIAS

Rotation matrix for the rotation from the International Celestial Reference Frame 
(`ICRF`) and the Mean Dynamical Equator and Equinox at J2000.0 (`MEME2000`).

### References
- Hilton, James L., and Catherine Y. Hohenkerk. -- Rotation matrix from the mean 
    dynamical equator and equinox at J2000. 0 to the ICRS. -- Astronomy & Astrophysics 
    413.2 (2004): 765-770. DOI: [10.1051/0004-6361:20031552](https://www.aanda.org/articles/aa/pdf/2004/02/aa3851.pdf)
- [SOFA docs](https://www.iausofa.org/2021_0512_C/sofa/sofa_pn_c.pdf)
"""
ROTMAT_ICRF2J2000_BIAS

# ------------------------------------------------------------------------------
#                               MEME2000 <-> ECLIPJ2000
# ------------------------------------------------------------------------------

function Rotation(
    origin::MeanEquatorMeanEquinoxJ2000, target::EclipticEquinoxJ2000, e::Epoch)
    return Rotation(origin, target, Orient.ROTMAT_J20002ECLIPJ2000)
end

function Rotation(
    origin::EclipticEquinoxJ2000, target::MeanEquatorMeanEquinoxJ2000, e::Epoch)
    return Rotation(origin, target, Orient.ROTMAT_J20002ECLIPJ2000')
end

connect!(MEME2000, ECLIPJ2000)

const ROTMAT_J20002ECLIPJ2000 = angle_to_dcm(
    Orient.orient_obliquity(Orient.iau2006a, 0.0), :X
)

"""
    ROTMAT_J20002ECLIPJ2000

Rotation matrix for the rotation from the Mean Dynamical Equator of J2000 
[`MEME2000`](@ref) to the Mean Ecliptic Equinox. 
This corresponds to the transformation `J2000 -> ECLIPJ2000` in the SPICE 
toolkit.
"""
ROTMAT_J20002ECLIPJ2000

# ------------------------------------------------------------------------------
#                               ICRF <-> MEMEMOD
# ------------------------------------------------------------------------------

function Rotation(
    origin::InternationalCelestialReferenceFrame, target::MeanEquatorMeanEquinoxOfDate, 
    e::Epoch; model::Orient.IAU2006Model=iau2006b)

    t = j2000c(convert(TT, e))
    γ, ϕ, ψ, ε = Orient.fw_angles(model, t)
    R = Orient.fw_matrix(γ, ϕ, ψ, ε)
    return Rotation(origin, target, R)
end

function Rotation(
    origin::MeanEquatorMeanEquinoxOfDate, 
    target::InternationalCelestialReferenceFrame, e::Epoch; 
    model::Orient.IAU2006Model=iau2006b)
    inv(Rotation(target, origin, e; model=model))
end

connect!(ICRF, MEMEMOD)