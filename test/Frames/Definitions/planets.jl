# load constants 
kclear()
constants = Basic.load(TPC(path(KERNELS[:IAU]))) 
furnsh(path(KERNELS[:IAU]))

# create dummy frame system 
FRAMES = FrameSystem{3, Float64}()

# insert ICRF
@axes ICRF 1 InternationalCelestialReferenceFrame
@axes IAU_TEST 2

@testset "Body-Centered Rotating (TOD)" verbose=true begin

    NAIFIds = []
    for k in keys(constants)
        (k>10 && k<1000 && haskey(constants[k], :pole_ra)) ? push!(NAIFIds, k) : ()
    end

    @testset "Body with NAIFId $i" for i in NAIFIds 

        FRAMES = FrameSystem{3, Float64}()
        add_axes_inertial!(FRAMES, ICRF)
        add_axes_bcrtod!(FRAMES, constants, i, "test", IAU_TEST, ICRF)

        for _ in 1:25
            ep = rand(0.0:1e6)

            R = sxform("J2000", "IAU_$(bodc2n(i))", ep)
            R_ = rotation6(FRAMES, ICRF, IAU_TEST, ep)

            v = rand(3)
            v /= norm(v)

            # It has a precision of about 36 marcsec
            @test acosd(max(-1, min(1, dot(R[1:3, 1:3]*v, R_[1]*v)))) ≤ 1e-5
            # @test acosd(max(-1, min(1, dot(R[4:6, 1:3]*v, R_[2]*v)))) ≤ 1e-5
            
        end
    end;

end;


# Test against manual computation 
p = Basic.Orient.PlanetsPrecessionNutation(601, constants)
f1, f2, f3 = Basic.Orient.orient_planets_angles(p, "mimas")

function mimas(sec)
    D = sec/Tempo.DAY2SEC
    T = D/Tempo.CENTURY2DAY
    S3 = 177.40 - 36505.5*T
    S5 = 316.45 + 506.2*T

    a = 40.66 - 0.036*T + 13.56*sind(S3)
    d = 83.52 - 0.004*T - 1.53*cosd(S3)
    w = 333.46 + 381.99455*D  - 13.48*sind(S3) - 44.85*sind(S5)

    return deg2rad.([a, d, w])
end

@testset "Euler angles and derivatives" begin
    
    for _ in 1:50
        e = rand(0.0:1e8)

        # angle 
        angles = [f1(e)...]
        angles_true = mimas(e)
        @test angles ≈ angles_true

        # 1st derivative 
        δangles = [f2(e)[4:end]...]
        δangles_true = ForwardDiff.derivative(mimas, e)
        @test δangles ≈ δangles_true

        # 2nd derivative 
        δ²angles = [f3(e)[7:end]...]
        δ²angles_true = ForwardDiff.derivative(t-> ForwardDiff.derivative(mimas, t), e)
        @test δ²angles ≈ δ²angles_true
    end

end