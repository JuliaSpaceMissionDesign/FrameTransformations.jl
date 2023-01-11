using BenchmarkTools
using ReferenceFrameRotations
using StaticArrays
using ForwardDiff
using Basic 
using Test

import Basic.Tempo: DJ2000


eph = CalcephProvider(["/home/michele/spice/kernels/spk/de440.bsp"])

# Frame with ephemeris loading 
FRAMES = FrameSystem{4, Float64}(eph);

icrf2meme = angle_to_dcm(pi/7, :X)
meme2eclip = angle_to_dcm(pi/4, pi/3, :ZY)
fun_itrf(t) = angle_to_dcm(t, :Z)

# Register axes! 
@axes ICRF 1 InternationalCelestialReferenceFrame
@axes MEME2000 2 MeanEarthMeanEquinoxJ2000
@axes ECLIPJ2000 3 EclipticEquinoxJ2000
@axes IAU_EARTH 4 

add_axes_inertial!(FRAMES, ICRF)
add_axes_inertial!(FRAMES, MEME2000; parent=ICRF, dcm=icrf2meme)
add_axes_inertial!(FRAMES, ECLIPJ2000; parent=MEME2000, dcm=meme2eclip)
add_axes_rotating!(FRAMES, IAU_EARTH, ICRF, fun_itrf)

# Register points! 
@point SSB 0 SolarSystemBarycenter
@point EMB 3 EarthMoonBarycenter 
@point Sun 10
@point Earth 399
@point VB 2 VenusBarycenter
@point Venus 299 Venus 

add_point_root!(FRAMES, SSB, ICRF)
add_point_ephemeris!(FRAMES, EMB)
add_point_ephemeris!(FRAMES, Earth)
add_point_ephemeris!(FRAMES, Sun)
add_point_ephemeris!(FRAMES, VB)
add_point_ephemeris!(FRAMES, Venus)

@benchmark vector3($FRAMES, $Venus, $SSB, $ICRF, 0.)
@benchmark vector6($FRAMES, $Venus, $EMB, $ICRF, 0.)
@benchmark rotation3($FRAMES, $IAU_EARTH, $ECLIPJ2000, 0.)


@testset "Point Transformations" verbose=true begin 
    fcns = (vector3, vector6, vector9, vector12)

    @testset "Identity Translations" verbose=true begin 
        for (i, fcn) in enumerate(fcns)
            stv = fcn(FRAMES, Sun, Sun, ECLIPJ2000, 0.)
            @test stv == zeros(3*i)
        end
    end

    jd0 = ephem_timespan(eph)[1] + 100
    t = 0.

    @testset "Ephemeris Translation" verbose=true begin 
        for (i, fcn) in enumerate(fcns)
            stv = fcn(FRAMES, Venus, SSB, ICRF, t)

            y = zeros(3i)
            ephem_compute_order!(y, eph, DJ2000, t, 0, 299, i-1)
            
            @test stv ≈ y atol=1e-10

        end
    end

    @testset "Rotated Translation" verbose=true begin 
        for (i, fcn) in enumerate(fcns)
            stv = fcn(FRAMES, Venus, SSB, ECLIPJ2000, t)

            y = zeros(3i)
            ephem_compute_order!(y, eph, DJ2000, t, 0, 299, i-1)
            
            R = Rotation(meme2eclip*icrf2meme, [DCM(0.0I) for i = 2:i]...)
            y = R*y 
            @test stv ≈ y atol=1e-10

        end
    end
end;


@testset "Rotation Transformations" verbose=true begin 
    fcns = (rotation3, rotation6, rotation9, rotation12)

    @testset "Identity Rotation" verbose=true begin 
        for (i, fcn) in enumerate(fcns)
            R = fcn(FRAMES, ICRF, ICRF, 0.)

            @test R[1] == DCM(1.0I)
            for k = 2:i 
                @test R[k] ≈ DCM(0.0I) atol=1e-11
            end
        end
    end

    @testset "1 Step Inertial Rotation" verbose=true begin
        for (i, fcn) in enumerate(fcns)    
            R1 = fcn(FRAMES, ICRF, MEME2000, 0.)
            R2 = fcn(FRAMES, MEME2000, ECLIPJ2000, 0.)

            for (j, R) in enumerate([R1, R2])
                Rₑ = j == 1 ? icrf2meme : meme2eclip
                @test R[1] ≈ Rₑ atol=1e-11
                for k = 2:i 
                    @test R[k] == DCM(0.0I) 
                end
            end
        end
    end;

    @testset "2 Step Inertial Rotation" verbose=true begin
        for (i, fcn) in enumerate(fcns)    
            R = fcn(FRAMES, ICRF, ECLIPJ2000, 0.)

            @test R[1] ≈ meme2eclip*icrf2meme atol=1e-11
            for k = 2:i 
                @test R[k] == DCM(0.0I)
            end

        end
    end;

    @testset "Rotating Rotation" verbose=true begin
        for (i, fcn) in enumerate(fcns)    
            t = π/3
            R = fcn(FRAMES, ECLIPJ2000, IAU_EARTH, t)

            @test R[1] ≈ fun_itrf(t)*icrf2meme'*meme2eclip' atol=1e-11

            if i >= 2  
                s, c = sincos(t)
                δ = DCM((-s, -c, 0, c, -s, 0, 0, 0, 0))

                @test R[2] ≈ δ*icrf2meme'*meme2eclip' atol=1e-11
                
                if i >= 3
                    δ² = DCM((-c, s, 0, -s, -c, 0, 0, 0, 0))
                    @test R[3] ≈ δ²*icrf2meme'*meme2eclip' atol=1e-11
                
                    if i >= 4 
                        δ³ = DCM((s, c, 0, -c, s, 0, 0, 0, 0))
                        @test R[4] ≈ δ³*icrf2meme'*meme2eclip' atol=1e-11
                
                    end
                end

            end

        end
    end;
end;
