import ReferenceFrameRotations: angle_to_dcm
import Basic.Utils: angle_to_δdcm, angle_to_δ²dcm, angle_to_δ³dcm
import Basic.Utils: D¹, D², D³

getdcm(t, seq, fx) = angle_to_dcm(fx(t), seq)
getdcm(t, seq, fx, fy) = angle_to_dcm(fx(t), fy(t), seq)
getdcm(t, seq, fx, fy, fz) = angle_to_dcm(fx(t), fy(t), fz(t), seq)

function getv(f, t)
    [f(t), D¹(f, t), D²(f, t), D³(f, t)]
end

fx = x->3x^2*cos(x)-sin(x)
fy = y->π/7*cos(y)
fz = z->π/3*sin(z)

@testset "DCM" verbose=true begin 

    τ = rand()
    x = getv(fx, τ)
    y = getv(fy, τ)
    z = getv(fz, τ)

    # 1 Rotation
    for seq in (:X, :Y, :Z)
        for (D, fun) in zip((D¹, D², D³), (angle_to_δdcm, angle_to_δ²dcm, angle_to_δ³dcm))
            d = D(t->getdcm(t, seq, fx), τ) 
            E = fun(x, seq)
            @test E ≈ d atol=1e-10
        end
    end

    # 2 Rotations 
    for seq in (:YX, :XY, :ZX, :XZ, :YZ, :ZY)
        for (D, fun) in zip((D¹, D²), (angle_to_δdcm, angle_to_δ²dcm))
            d = D(t->getdcm(t, seq, fx, fy), τ) 
            E = fun(x, y, seq)
            @test E ≈ d atol=1e-10
        end
    end

    # 3 Rotations
    for seq in (:ZYX, :XYX, :XYZ, :XZX, :XZY, :YXY, 
                 :YXZ, :YZX, :YZY, :ZXY, :ZXZ, :ZYZ)

        for (D, fun) in zip((D¹, D²), (angle_to_δdcm, angle_to_δ²dcm))
            d = D(t->getdcm(t, seq, fx, fy, fz), τ) 
            E = fun(x, y, z, seq)
            @test E ≈ d atol=1e-10
        end
    end


end;
