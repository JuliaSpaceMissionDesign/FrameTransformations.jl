export orient_obliquity


"""
    orient_obliquity(m::IAUModel, t::Number)

Compute the mean obliquity of the ecliptic at epoch, in radians, at time `t` expressed 
in `TT` Julian centuries since [`J2000`](@ref). 

!!! note 
    This function is implemented only for `IAU1980` and `IAU2006` models. IAU 2000 Models 
    implement proper precession-rate corrections to the IAU1980 mean obliquity. 

### References 
- ERFA [obl80](https://github.com/liberfa/erfa/blob/master/src/obl80.c) and 
  [obl06](https://github.com/liberfa/erfa/blob/master/src/obl06.c) functions.
"""
function orient_obliquity(::IAU2006Model, t::Number)
    return @evalpoly(t,
        84381.406, -46.836769, -0.0001831, 
        0.00200340, -0.000000576, -0.0000000434) |> arcsec2rad
end

function orient_obliquity(::IAU1980Model, t::Number)
    return @evalpoly(t,
        84381.448, -46.8150, -0.00059, 0.001813
    ) |> arcsec2rad
end

