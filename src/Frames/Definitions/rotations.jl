function Rotation(origin::ICRF, target::MEME2000, e::Epoch)
    return Rotation(origin, target, ICRF2J2000_BIAS)
end