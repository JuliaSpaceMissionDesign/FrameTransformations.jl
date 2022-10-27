function Rotation(
    origin::InternationalCelestialReferenceFrame,
    target::MeanEclipticOfDate, 
    e::Epoch)
    m = orient_icrs2ecleod(DJ2000, j2000(convert(TT, e)))
    Rotation(origin, target, m)
end

function Rotation(
    origin::MeanEclipticOfDate, 
    target::InternationalCelestialReferenceFrame, 
    e::Epoch)
    inv(Rotation(target, origin, e))
end

function Rotation(
    origin::InternationalCelestialReferenceFrame, 
    target::MeanEclipticJ2000, 
    e::Epoch)
    m = orient_icrs2ecleod(DJ2000, 0.0)
    Rotation(origin, target, m)
end

function Rotation(
    origin::MeanEclipticJ2000, 
    target::InternationalCelestialReferenceFrame, 
    e::Epoch)
    inv(Rotation(target, origin, e))
end