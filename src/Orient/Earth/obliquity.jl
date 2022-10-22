import Basic.Utils: sec2rad

"""
    obliquity(::IAU2006, jd1tt::N, jd2tt::N2) where {N<:Number, N2<:Number}

Mean obliquity of the ecliptic, IAU 2006 precession model.

### Inputs
- `jd1tt`, `jd2tt` --  TT as a 2-part Julian Date 

### Output
Obliquity of the ecliptic at epoch -- `rad`

### References 
- From [erfa](https://github.com/liberfa/erfa/blob/master/src/obl06.c) library,
"""
function obliquity(::IAU2006, jd1tt::N, jd2tt::N2) where {N<:Number, N2<:Number}
    #Interval between fundamental date J2000.0 and given date (JC).
    t = ((jd1tt - DJ2000) + jd1tt) / 36525.0;
    return @evalpoly(t,
        84381.406, -46.836769, -0.0001831, 
        0.00200340, -0.000000576, -0.0000000434) |> sec2rad
end

