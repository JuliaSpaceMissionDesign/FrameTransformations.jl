

function alg1(t2060::N, tt::N, long::N) where N
    wtt = 0.017202786*tt
    s1, c1 = sincos(wtt)
    s2 = 2.0*s1*c1
    c2 = (c1+s1)*(c1-s1)

    ra = -1.38880 + 1.72027920e-2*tt + 3.199e-2*s1 - 2.65e-3*c1 + 4.050e-2*s2 + 1.525e-2*c2
    ra = mod(ra, 2π)
    dec = 6.57e-3 + 7.347e-2*s1 - 3.9919e-1*c1 + 7.3e-4*s2 - 6.60e-3*c2
    ha = 1.75283 + 6.3003881*t2060 + long - ra
    ha = mod(ha + π, 2π) - π
    return ra, dec, ha

end



const pi2 = 2*pi
const pim = pi/2

function sunpos(dt::DateTime, deltaT::Float64,
                longitude::Float64, latitude::Float64, height::Float64=0.0,
                pressure::Float64=1.0, temperature::Float64=20.0; flag='l', alg=alg5)
    t2060,tt = date2060(dt, deltaT)
    rightAscension,declination,hourAngle = alg(t2060, tt, deltaT, longitude)
    
    if flag=='s'
        return rightAscension,declination,hourAngle
    end

    zenith,azimuth = zenithazimuth(latitude, declination, hourAngle, height,
                                   pressure, temperature)

    return flag=='l' ? (zenith,azimuth,rightAscension,declination,hourAngle) : (zenith,azimuth)
end # function sunpos


# with original time spec
function sunpos(ut::Float64, day::Int32, month::Int32, year::Int32, deltaT::Float64,
                longitude::Float64, latitude::Float64, height::Float64=0.0,
                pressure::Float64=1.0, temperature::Float64=20.0; flag='l', alg=alg5)
    t2060,tt = date2060(ut, day, month, year, deltaT)
    rightAscension,declination,hourAngle = alg(t2060, tt, deltaT, longitude)

    if flag=='s'
        return rightAscension,declination,hourAngle
    end

    zenith,azimuth = zenithazimuth(latitude, declination, hourAngle, height,
                                   pressure, temperature)

    return flag=='l' ? (zenith,azimuth,rightAscension,declination,hourAngle) : (zenith,azimuth)
end # function sunpos_alg1

const dt2060 = DateTime(2060,1,1)

deltaTlin(t2060::Float64) = 96.4 + 0.00158*t2060

# Time computation aka julia
function date2060(dt::DateTime, deltaT::Float64)
    t2060 = Int64(Dates.value(dt.instant-dt2060.instant))/86400000.0
    if isnan(deltaT)
        deltaT = deltaTlin(t2060)
    end

    return t2060,t2060 + 1.1574e-5*deltaT
end # function date2060

# Transcript of the original time scale computation:
function date2060(ut::Float64, day::Int32, month::Int32, year::Int32, deltaT::Float64)
    # number of days starting from the beginning of the year 2060:
    if month<=2
        mt = month+12
        yt = year-1
    else
        mt = month
        yt = year
    end

    t2060 = (itrunc(365.25*(yt-2000)) + itrunc(30.6001*(mt+1)) - itrunc(0.01*(yt)) + day) +
             0.0416667*ut - 21958.0
    if isnan(deltaT)
        deltaT = deltaTlin(t2060)
    end
    tt = t2060 + 1.1574e-5*deltaT

    return t2060,t2060 + 1.1574e-5*deltaT
end # function date2060

function zenithazimuth(latitude::Float64, declination::Float64, hourAngle::Float64,
                       height::Float64, pressure::Float64, temperature::Float64)
    sp = sin(latitude)
    cp = sqrt((1-sp*sp))
    sd = sin(declination)
    cd = sqrt(1-sd*sd)
    sH = sin(hourAngle)
    cH = cos(hourAngle)
    se0 = sp*sd + cp*cd*cH
    # elevation including parallax correction with height in km
    ep = asin(se0) - (6371.0 + height)/149597871*sqrt(1.0-se0*se0)
    azimuth = atan(sH, cH*sp - sd*cp/cd)

    zenith  = pim - ep
    zenith -= ep>0.0 ?
             (0.08422*pressure)/((273.0+temperature)*tan(ep + 0.003138/(ep + 0.08919))) : 0.0
    
    return zenith,azimuth
end # function zenithazimuth




function alg1(t2060::Float64, tt::Float64, deltaT::Float64, longitude::Float64)
    wtt = 0.017202786*tt

    s1 = sin(wtt)
    c1 = cos(wtt)
    s2 = 2.0*s1*c1
    c2 = (c1+s1)*(c1-s1)

    rightAscension = -1.38880 + 1.72027920e-2*tt + 3.199e-2*s1 - 2.65e-3*c1 +
                      4.050e-2*s2 + 1.525e-2*c2;
    rightAscension = mod(rightAscension, pi2)

    declination = 6.57e-3 + 7.347e-2*s1 - 3.9919e-1*c1 + 7.3e-4*s2 - 6.60e-3*c2

    hourAngle = 1.75283 + 6.3003881*t2060 + longitude - rightAscension
    hourAngle = mod(hourAngle + pi, pi2) - pi

    return rightAscension,declination,hourAngle
end # function alg1
                    
function sunpos_alg1(dt::DateTime, deltaT::Float64,
                     longitude::Float64, latitude::Float64, height::Float64=0.0,
                     pressure::Float64=1.0, temperature::Float64=20.0; flag='l')
    t2060,tt = date2060(dt, deltaT)
    rightAscension,declination,hourAngle = alg1(t2060, tt, deltaT, longitude)
    
    if flag=='s'
        return rightAscension,declination,hourAngle
    end

    zenith,azimuth = zenithazimuth(latitude, declination, hourAngle, height,
                                   pressure, temperature)

    return flag=='l' ? (zenith,azimuth,rightAscension,declination,hourAngle) : (zenith,azimuth)
end  # function sunpos_alg1

# with original time spec
function sunpos_alg1(ut::Float64, day::Int32, month::Int32, year::Int32, deltaT::Float64,
                     longitude::Float64, latitude::Float64, height::Float64=0.0,
                     pressure::Float64=1.0, temperature::Float64=20.0; flag='l')
    t2060,tt = date2060(ut, day, month, year, deltaT)
    rightAscension,declination,hourAngle = alg1(t2060, tt, deltaT, longitude)

    if flag=='s'
        return rightAscension,declination,hourAngle
    end

    zenith,azimuth = zenithazimuth(latitude, declination, hourAngle, height,
                                   pressure, temperature)

    return flag=='l' ? (zenith,azimuth,rightAscension,declination,hourAngle) : (zenith,azimuth)
end # function sunpos_alg1

function alg2(t2060::Float64, tt::Float64, deltaT::Float64, longitude::Float64)
    wtt = 0.017202786*tt

    s1 = sin(wtt)
    c1 = cos(wtt)
    s2 = 2.0*s1*c1
    c2 = (c1+s1)*(c1-s1)
    s3 = s2*c1 + c2*s1
    c3 = c2*c1 - s2*s1
    s4 = 2.0*s2*c2
    c4 = (c2+s2)*(c2-s2)

    rightAscension = -1.38880 + 1.72027920e-2*tt + 3.199e-2*s1 - 2.65e-3*c1 +
                      4.050e-2*s2 + 1.525e-2*c2 + 1.33e-3*s3 + 3.8e-4*c3 +
                      7.3e-4*s4 + 6.2e-4*c4
    rightAscension = mod(rightAscension, pi2)

    declination = 6.57e-3 + 7.347e-2*s1 - 3.9919e-1*c1 + 7.3e-4*s2 - 6.60e-3*c2 +
                  1.50e-3*s3 - 2.58e-3*c3 + 6e-5*s4 - 1.3e-4*c4

    hourAngle = 1.75283 + 6.3003881*t2060 + longitude - rightAscension
    hourAngle = mod(hourAngle + pi, pi2) - pi

    return rightAscension,declination,hourAngle
end # function alg2
                    
function sunpos_alg2(dt::DateTime, deltaT::Float64,
                     longitude::Float64, latitude::Float64, height::Float64=0.0,
                     pressure::Float64=1.0, temperature::Float64=20.0; flag='l')
    t2060,tt = date2060(dt, deltaT)
    rightAscension,declination,hourAngle = alg2(t2060, tt, deltaT, longitude)
    
    if flag=='s'
        return rightAscension,declination,hourAngle
    end

    zenith,azimuth = zenithazimuth(latitude, declination, hourAngle, height,
                                   pressure, temperature)

    return flag=='l' ? (zenith,azimuth,rightAscension,declination,hourAngle) : (zenith,azimuth)
end  # function sunpos_alg2

# with original time spec
function sunpos_alg2(ut::Float64, day::Int32, month::Int32, year::Int32, deltaT::Float64,
                     longitude::Float64, latitude::Float64, height::Float64=0.0,
                     pressure::Float64=1.0, temperature::Float64=20.0; flag='l')
    t2060,tt = date2060(ut, day, month, year, deltaT)
    rightAscension,declination,hourAngle = alg2(t2060, tt, deltaT, longitude)

    if flag=='s'
        return rightAscension,declination,hourAngle
    end

    zenith,azimuth = zenithazimuth(latitude, declination, hourAngle, height,
                                   pressure, temperature)

    return flag=='l' ? (zenith,azimuth,rightAscension,declination,hourAngle) : (zenith,azimuth)
end # function sunpos_alg2

function alg3(t2060::Float64, tt::Float64, deltaT::Float64, longitude::Float64)
    wtt = 0.017202786*tt

    lambda = -1.388803 + 1.720279216e-2*tt + 3.3366e-2*sin(wtt - 0.06172) +
              3.53e-4*sin(2.0*wtt - 0.1163)

    epsi = 4.089567e-1 - 6.19e-9*tt

    sl = sin(lambda)
    cl = cos(lambda)
    se = sin(epsi)
    ce = sqrt(1-se*se)

    rightAscension = atan(sl*ce, cl)
    if rightAscension<0.0 
        rightAscension += pi2
    end

    declination = asin(sl*se)

    hourAngle = 1.7528311 + 6.300388099*t2060 + longitude - rightAscension
    hourAngle = mod(hourAngle + pi, pi2) - pi;

    return rightAscension,declination,hourAngle
end # function alg3
     
function sunpos_alg3(dt::DateTime, deltaT::Float64,
                     longitude::Float64, latitude::Float64, height::Float64=0.0,
                     pressure::Float64=1.0, temperature::Float64=20.0; flag='l')
    t2060,tt = date2060(dt, deltaT)
    rightAscension,declination,hourAngle = alg3(t2060, tt, deltaT, longitude)
    
    if flag=='s'
        return rightAscension,declination,hourAngle
    end

    zenith,azimuth = zenithazimuth(latitude, declination, hourAngle, height,
                                   pressure, temperature)

    return flag=='l' ? (zenith,azimuth,rightAscension,declination,hourAngle) : (zenith,azimuth)
end  # function sunpos_alg3

# with original time spec
function sunpos_alg3(ut::Float64, day::Int32, month::Int32, year::Int32, deltaT::Float64,
                     longitude::Float64, latitude::Float64, height::Float64=0.0,
                     pressure::Float64=1.0, temperature::Float64=20.0; flag='l')
    t2060,tt = date2060(ut, day, month, year, deltaT)
    rightAscension,declination,hourAngle = alg3(t2060, tt, deltaT, longitude)

    if flag=='s'
        return rightAscension,declination,hourAngle
    end

    zenith,azimuth = zenithazimuth(latitude, declination, hourAngle, height,
                                   pressure, temperature)

    return flag=='l' ? (zenith,azimuth,rightAscension,declination,hourAngle) : (zenith,azimuth)
end # function sunpos_alg3

function alg4(t2060::Float64, tt::Float64, deltaT::Float64, longitude::Float64)
    wtt = 0.017202786*tt

    L = 1.752790 + 1.720279216e-2*tt + 3.3366e-2*sin(wtt - 0.06172) +
        3.53e-4*sin(2.0*wtt - 0.1163)

    nu = 9.282e-4*tt - 0.8
    dlam = 8.34e-5*sin(nu)
    lambda = L + pi + dlam

    epsi = 4.089567e-1 - 6.19e-9*tt + 4.46e-5*cos(nu)

    sl = sin(lambda)
    cl = cos(lambda)
    se = sin(epsi)
    ce = sqrt(1-se*se)

    rightAscension = atan(sl*ce, cl);
    if rRightAscension<0.0
        rightAscension += pi2
    end

    declination = asin(sl*se)

    hourAngle = 1.7528311 + 6.300388099*t2060 + longitude - rightAscension + 0.92*dlam;
    hourAngle = mod(hourAngle + pi, pi2) - pi;

    return rightAscension,declination,hourAngle
end # function alg4
     
function sunpos_alg4(dt::DateTime, deltaT::Float64,
                     longitude::Float64, latitude::Float64, height::Float64=0.0,
                     pressure::Float64=1.0, temperature::Float64=20.0; flag='l')
    t2060,tt = date2060(dt, deltaT)
    rightAscension,declination,hourAngle = alg4(t2060, tt, deltaT, longitude)
    
    if flag=='s'
        return rightAscension,declination,hourAngle
    end

    zenith,azimuth = zenithazimuth(latitude, declination, hourAngle, height,
                                   pressure, temperature)

    return flag=='l' ? (zenith,azimuth,rightAscension,declination,hourAngle) : (zenith,azimuth)
end  # function sunpos_alg4

# with original time spec
function sunpos_alg4(ut::Float64, day::Int32, month::Int32, year::Int32, deltaT::Float64,
                     longitude::Float64, latitude::Float64,  height::Float64=0.0,
                     pressure::Float64=1.0, temperature::Float64=20.0; flag='l')
    t2060,tt = date2060(ut, day, month, year, deltaT)
    rightAscension,declination,hourAngle = alg4(t2060, tt, deltaT, longitude)

    if flag=='s'
        return rightAscension,declination,hourAngle
    end

    zenith,azimuth = zenithazimuth(latitude, declination, hourAngle, height,
                                   pressure, temperature)

    return flag=='l' ? (zenith,azimuth,rightAscension,declination,hourAngle) : (zenith,azimuth)
end # function sunpos_alg4

function alg5(t2060::Float64, tt::Float64, deltaT::Float64, longitude::Float64)
    wtt = 0.0172019715*tt

    s1 = sin(wtt)
    c1 = cos(wtt)
    s2 = 2.0*s1*c1
    c2 = (c1+s1)*(c1-s1)
    s3 = s2*c1 + c2*s1
    c3 = c2*c1 - s2*s1

    L = 1.7527901 + 1.7202792159e-2*tt + 3.33024e-2*s1 - 2.0582e-3*c1 +
        3.512e-4*s2 - 4.07e-5*c2 + 5.2e-6*s3 - 9e-7*c3 -
        8.23e-5*s1*sin(2.92e-5*tt) + 1.27e-5*sin(1.49e-3*tt - 2.337) +
        1.21e-5*sin(4.31e-3*tt + 3.065) + 2.33e-5*sin(1.076e-2*tt - 1.533) +
        3.49e-5*sin(1.575e-2*tt - 2.358) + 2.67e-5*sin(2.152e-2*tt + 0.074) +
        1.28e-5*sin(3.152e-2*tt + 1.547) + 3.14e-5*sin(2.1277e-1*tt - 0.488)

    nu = 9.282e-4*tt - 0.8
    dlam = 8.34e-5*sin(nu)
    lambda = L + pi + dlam

    epsi = 4.089567e-1 - 6.19e-9*t2060 + 4.46e-5*cos(nu)

    sl = sin(lambda)
    cl = cos(lambda)
    se = sin(epsi)
    ce = sqrt(1-se*se)

    rightAscension = atan(sl*ce, cl)
    if rightAscension<0.0
        rightAscension += pi2
    end

    declination = asin(sl*se)

    hourAngle = 1.7528311 + 6.300388099*t2060 + longitude - rightAscension + 0.92*dlam;
    hourAngle = mod(hourAngle + pi, pi2) - pi;

    return rightAscension,declination,hourAngle
end # function alg5
     
function sunpos_alg5(dt::DateTime, deltaT::Float64,
                     longitude::Float64, latitude::Float64, height::Float64=0.0,
                     pressure::Float64=1.0, temperature::Float64=20.0; flag='l')
    t2060,tt = date2060(dt, deltaT)
    rightAscension,declination,hourAngle = alg5(t2060, tt, deltaT, longitude)
    
    if flag=='s'
        return rightAscension,declination,hourAngle
    end

    zenith,azimuth = zenithazimuth(latitude, declination, hourAngle, height,
                                   pressure, temperature)

    return flag=='l' ? (zenith,azimuth,rightAscension,declination,hourAngle) : (zenith,azimuth)
end  # function sunpos_alg5

# with original time spec
function sunpos_alg5(ut::Float64, day::Int32, month::Int32, year::Int32, deltaT::Float64,
                     longitude::Float64, latitude::Float64, height::Float64=0.0,
                     pressure::Float64=1.0, temperature::Float64=20.0; flag='l')
    t2060,tt = date2060(ut, day, month, year, deltaT)
    rightAscension,declination,hourAngle = alg5(t2060, tt, deltaT, longitude)

    if flag=='s'
        return rightAscension,declination,hourAngle
    end

    zenith,azimuth = zenithazimuth(latitude, declination, hourAngle, height,
                                   pressure, temperature)

    return flag=='l' ? (zenith,azimuth,rightAscension,declination,hourAngle) : (zenith,azimuth)
end # function sunpos_alg4
