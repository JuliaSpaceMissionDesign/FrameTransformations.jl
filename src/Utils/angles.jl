export sec2rad, rad2secW

const ARCSECONDS_PER_DEGREE = 3600

@inline sec2rad(x) = deg2rad(x / ARCSECONDS_PER_DEGREE)

@inline rad2sec(x) = rad2deg(x) * ARCSECONDS_PER_DEGREE