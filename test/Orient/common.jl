
function _test_matrix(e, g, fun, tol=1e-8)
    v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

    ft = fun(j2000s(e))
    err = g*ft' - I(3)
    @test norm(err) ≤ tol

    for i in 1:3 
        v = zeros(3)
        v[i] = 1.0
        @test v2as(g * v, ft * v) ≤ tol
        v[i] = -1.0
        @test v2as(g * v, ft * v) ≤ tol
    end
end

@testset "ICRF-Based (GODOT)" verbose=true begin

    @testset "MOD" begin
        # Performed against GODOT v1.4.0

        e = Epoch("2007-01-10T17:48:30.204000 TT")
        g = DCM(
            9.99998533e-01, -1.57114117e-03, -6.82595764e-04,
            1.57114120e-03,  9.99998766e-01, -4.94314071e-07,
            6.82595698e-04, -5.78140984e-07,  9.99999767e-01)' 
        _test_matrix(e, g, Orient.orient_rot3_icrf_to_mod)

        e = Epoch("2016-01-09T13:59:54.204000 TT")
        g = DCM(
            9.99992370e-01, -3.58284487e-03, -1.55666798e-03,
            3.58284495e-03,  9.99993582e-01, -2.73547078e-06,
            1.55666779e-03, -2.84185011e-06,  9.99998788e-01)'
        _test_matrix(e, g, Orient.orient_rot3_icrf_to_mod)

        e = Epoch("2022-01-01T12:00:00.0 TT")
        g = DCM(
            9.99985612e-01, -4.91996438e-03, -2.13759287e-03, 
            4.91996451e-03,  9.99987897e-01, -5.19786808e-06, 
            2.13759257e-03, -5.31908778e-06,  9.99997715e-01)'
        _test_matrix(e, g, Orient.orient_rot3_icrf_to_mod)

        e = Epoch("2029-07-08T08:17:00.204000 TT")
        g = DCM(
            9.99974104e-01, -6.60056979e-03, -2.86769780e-03, 
            6.60056999e-03,  9.99978216e-01, -9.39453283e-06, 
            2.86769734e-03, -9.53415048e-06,  9.99995888e-01)'
        _test_matrix(e, g, Orient.orient_rot3_icrf_to_mod)

    end

    @testset "TOD" begin
        # Performed against GODOT v1.4.0

        e = Epoch("2026-11-10T22:14:26.196000 TT")
        g = DCM(
            9.99978303e-01, -6.04182243e-03, -2.62493419e-03,
            6.04172827e-03,  9.99981748e-01, -4.38006829e-05,
            2.62515091e-03,  2.79405935e-05,  9.99996554e-01)'
        _test_matrix(e, g, Orient.orient_rot3_icrf_to_tod)

        e = Epoch("2015-01-19T17:05:52.332000 TT")
        g = DCM(
            9.99993168e-01, -3.39019667e-03, -1.47294334e-03,
            3.39026374e-03,  9.99994252e-01,  4.30407853e-05,
            1.47278896e-03, -4.80341577e-05,  9.99998914e-01)'
        _test_matrix(e, g, Orient.orient_rot3_icrf_to_tod)
        
        e = Epoch("2002-02-06T04:42:39.600000 TT")
        g = DCM(
            9.99999905e-01, -3.99160254e-04, -1.73409774e-04,
            3.99159426e-04,  9.99999920e-01, -4.80651508e-06,
            1.73411679e-04,  4.73729648e-06,  9.99999985e-01)'
        _test_matrix(e, g, Orient.orient_rot3_icrf_to_tod)
        
    end

    @testset "TEME" begin
        # Performed against GODOT v1.4.0
        tol = 1e-3

        e = Epoch("2025-01-28T05:38:10.536000 TT")
        g = DCM(
            9.99981304e-01, -5.60738226e-03, -2.43907997e-03,
            5.60727709e-03,  9.99984278e-01, -4.99554106e-05,
            2.43932174e-03,  3.62778794e-05,  9.99997024e-01)'
        err = g*Orient.orient_rot3_icrf_to_teme(j2000s(e))' - I(3)
        @test norm(err) ≤ tol

        e = Epoch("2029-08-22T06:02:21.552000 TT")
        g = DCM(
            9.99973784e-01, -6.62806523e-03, -2.91532267e-03,
            6.62805685e-03,  9.99978034e-01, -1.25351237e-05,
            2.91534172e-03, -6.78812935e-06,  9.99995750e-01)'
        err = g*Orient.orient_rot3_icrf_to_teme(j2000s(e))' - I(3)
        @test norm(err) ≤ tol
    end

end

@testset "IAU-76/FK5 Reduction (Vallado)" verbose=true begin
    # Based on Vallado (2013), pag. 230, example 3-15

    e = Epoch("2004-04-06T07:51:28.386009 UTC")
    e = convert(TT, e)

    r_itrf = SA[-1033.4793830, 7901.2952754, 6380.3565958]
    v_itrf = SA[-3.225636520, -2.872451450, 5.531924446]

    r_pef = SA[-1033.4750313, 7901.3055856, 6380.3445328]
    v_pef = SA[-3.225632747,  -2.872442511, 5.531931288]

    r_tod = SA[5094.5147804, 6127.3664612, 6380.3445328]
    v_tod = SA[-4.746088567, 0.786077222, 5.531931288]

    r_mod = SA[5094.0283745, 6127.8708164, 6380.2485164]
    v_mod = SA[-4.746263052, 0.786014045, 5.531790562]

    r_gcrf = SA[5102.508958, 6123.011401, 6378.136928]
    v_gcrf = SA[-4.74322016, 0.79053650, 5.53375528]

    @testset "ITRF to PEF" begin
        R_itrf2pef = Orient.orient_rot3_itrf_to_pef(j2000s(e))
        @test all(r_pef - R_itrf2pef * r_itrf .≤ 1e-4)  
        @show all(v_pef - R_itrf2pef * v_itrf .≤ 1e-6)  
    end

    @testset "PEF to TOD" begin
        R_pef2tod, R_pef2tod_d = Orient.orient_rot6_pef_to_tod(j2000s(e))
        @test all(r_tod - R_pef2tod * r_pef .≤ 1e-4)  

        GAST_ref = 312.8067654 # deg
        GAST = Orient.orient_gast(iau2000a, j2000s(e)/Tempo.CENTURY2SEC) |> rad2deg
        @test isapprox(GAST_ref, GAST, atol=1e-6)

        # Velocity correction 
        wxr_pef = SA[-0.57617229, -0.07536219, 0.0] 
        @test all(wxr_pef - R_pef2tod_d * r_pef .≤ 1e-6)

        v = R_pef2tod * (v_pef + R_pef2tod_d * r_pef)
        @test all(v_tod - v .≤ 1e-6)  
    end

    @testset "TOD to MOD" begin
        R_tod2mod = Orient.orient_rot3_tod_to_mod(j2000s(e))
        # TODO: km error here, due to iau2006 model used instead of 1980, update when 1980 avaliable
        @test all(r_mod - R_tod2mod * r_tod .≤ 1.0)  
        @test all(v_mod - R_tod2mod * v_tod .≤ 1e-3)  
    end

    @testset "ICRF/GCRF to TOD" begin
        @test all(r_tod - Orient.orient_rot3_icrf_to_tod(j2000s(e)) * r_gcrf .≤ 1e-3)
        @test all(v_tod - Orient.orient_rot3_icrf_to_tod(j2000s(e)) * v_gcrf .≤ 1e-6)  
    end

    @testset "ICRF/GCRF to MOD" begin
        @test all(r_mod - Orient.orient_rot3_icrf_to_mod(j2000s(e)) * r_gcrf .≤ 1e-3)
        @test all(v_mod - Orient.orient_rot3_icrf_to_mod(j2000s(e)) * v_gcrf .≤ 1e-6)  
    end

    @testset "ICRF/GCRF to MOD" begin
        @test all(r_mod - Orient.orient_rot3_icrf_to_mod(j2000s(e)) * r_gcrf .≤ 1e-3)
        @test all(v_mod - Orient.orient_rot3_icrf_to_mod(j2000s(e)) * v_gcrf .≤ 1e-6)  
    end
    
end
