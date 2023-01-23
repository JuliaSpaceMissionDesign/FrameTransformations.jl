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

            # @test isapprox(dot(R[1:3, 1:3] *v, R_[1]*v), 1.0, atol=1e-8)
            
        end
    end;

end;


