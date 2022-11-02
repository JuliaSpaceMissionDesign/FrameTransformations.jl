# ------------------------------------------------------------------------------
#                               ICRF <-> MEME2000
# ------------------------------------------------------------------------------

function Rotation(
    origin::InternationalCelestialReferenceFrame, 
    target::MeanEquatorMeanEquinoxJ2000, e::Epoch)
    return Rotation(origin, target, Orient.ICRF2J2000_BIAS)
end

function Rotation(
    origin::MeanEquatorMeanEquinoxJ2000, 
    target::InternationalCelestialReferenceFrame, e::Epoch)
    Rotation(origin, target, Orient.ICRF2J2000_BIAS')
end

connect!(ICRF, MEME2000)

# ------------------------------------------------------------------------------
#                               MEME2000 <-> ECLIPJ2000
# ------------------------------------------------------------------------------

function Rotation(
    origin::MeanEquatorMeanEquinoxJ2000, 
    target::EclipticEquinoxJ2000, e::Epoch)
    return Rotation(origin, target, Orient.EQ2ECL_J2000)
end

function Rotation(
    origin::EclipticEquinoxJ2000, 
    target::MeanEquatorMeanEquinoxJ2000, e::Epoch)
    return Rotation(origin, target, Orient.EQ2ECL_J2000')
end

connect!(MEME2000, ECLIPJ2000)

# ------------------------------------------------------------------------------
#                               ICRF <-> MEMEOD
# ------------------------------------------------------------------------------

function Rotation(
    origin::InternationalCelestialReferenceFrame, 
    target::MeanEquatorMeanEquinoxOfDate, e::Epoch; 
    model::Orient.IAU2006Model=iau2006b)

    t = j2000c(convert(TT, e))
    γ, ϕ, ψ, ε = Orient.fw_angles(model, t)
    PB = angle_to_dcm(γ, :Z) * angle_to_dcm(ϕ, -ψ, -ε, :XZX)
    E = angle_to_dcm(ε, :X)
    Rmat = E * PB

    return Rotation(origin, target, Rmat)
end

function Rotation(
    origin::MeanEquatorMeanEquinoxOfDate, 
    target::InternationalCelestialReferenceFrame, e::Epoch; 
    model::Orient.IAU2006Model=iau2006b)
    inv(Rotation(target, origin, e; model=model))
end

connect!(ICRF, MEMEOD)