using Basic.Tempo

function _random_datetime()
    ry = rand(1800:2100)
    rm = rand(1:12)
    ly = isleapyear(ry)
    if ly && rm == 2
        rd = rand(1:29)
    elseif rm == 2
        rd = rand(1:28)
    else
        rd = rand(1:Tempo.MTAB[rm])
    end
    rH = rand(0:23)
    rM = rand(0:59)
    rS = rand(0:59)
    rF = rand(0.0:0.000000001:1.0)

    return ry, rm, rd, rH, rM, rS, rF
end

function _random_datetime_isostr()
    ry, rm, rd, rH, rM, rS, rF = _random_datetime()
    rs = rS + rF
    return "$ry-$rm-$(rd)T$rH:$rM:$rs", ry, rm, rd, rH, rM, rS, rF
end

function _random_epoch()
    _, ry, rm, rd, rH, rM, rS, rF = _random_datetime_isostr()
    rs = rS + rF
    nscales = length(Tempo.TIMESCALES_ACRONYMS)
    rscale = rand(1:nscales)
    return "$ry-$rm-$(rd)T$rH:$rM:$rs $(String(Tempo.TIMESCALES_ACRONYMS[rscale]))",
    ry, rm, rd, rH, rM, rS,
    rF
end

include("convert.jl")
include("parse.jl")
include("offset.jl")
include("scales.jl")
include("datetime.jl")
include("origin.jl")
include("epoch.jl")
