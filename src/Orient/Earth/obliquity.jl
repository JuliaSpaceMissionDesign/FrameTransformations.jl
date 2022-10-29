export orient_obliquity

"""
    orient_obliquity(::M, j2000ttc::N) where {N<:Number, M<:IAU2006Model}

Mean obliquity of the ecliptic, IAU 2006 precession model.

### Inputs
- `t` --  `TT` centuries since J2000

### Output
Obliquity of the ecliptic at epoch -- `rad`

### References 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/obl06.c) library
"""
function orient_obliquity(::M, t::N) where {N<:Number, M<:IAU2006Model}
    return @evalpoly(t,
        84381.406, -46.836769, -0.0001831, 
        0.00200340, -0.000000576, -0.0000000434) |> arcsec2rad
end

const EQ2ECL_J2000 = angle_to_dcm(
    orient_obliquity(iau2006a, 0.0), :X
)