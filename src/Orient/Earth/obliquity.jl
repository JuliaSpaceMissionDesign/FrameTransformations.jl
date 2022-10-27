"""
    mean_obliquity(::M, j2000ttc::N) where {N<:Number, M<:IAU2006Model}

Mean obliquity of the ecliptic, IAU 2006 precession model.

### Inputs
- `j2000ttc` --  `TT` centuries since J2000

### Output
Obliquity of the ecliptic at epoch -- `rad`

### References 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/obl06.c) library
"""
function mean_obliquity(::M, j2000ttc::N) where {N<:Number, M<:IAU2006Model}
    return @evalpoly(j2000ttc,
        84381.406, -46.836769, -0.0001831, 
        0.00200340, -0.000000576, -0.0000000434) |> arcsec2rad
end

