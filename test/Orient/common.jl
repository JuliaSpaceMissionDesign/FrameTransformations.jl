@testset "MOD" begin
    # Performed against GODOT v1.4.0

    e = Epoch("2007-01-10T17:48:30.204000 TT")
    m0 = DCM(
        9.99998533e-01, -1.57114117e-03, -6.82595764e-04,
        1.57114120e-03,  9.99998766e-01, -4.94314071e-07,
        6.82595698e-04, -5.78140984e-07,  9.99999767e-01)' 

    @test all(m0 .≈ Orient.orient_rot3_icrf_to_mod(j2000s(e)))

    e = Epoch("2016-01-09T13:59:54.204000 TT")
    m0 = DCM(
        9.99992370e-01, -3.58284487e-03, -1.55666798e-03,
        3.58284495e-03,  9.99993582e-01, -2.73547078e-06,
        1.55666779e-03, -2.84185011e-06,  9.99998788e-01)'

    @test all(m0 .≈ Orient.orient_rot3_icrf_to_mod(j2000s(e)))


    e = Epoch("2022-01-01T12:00:00.0 TT")
    m0 = DCM(
        9.99985612e-01, -4.91996438e-03, -2.13759287e-03, 
        4.91996451e-03,  9.99987897e-01, -5.19786808e-06, 
        2.13759257e-03, -5.31908778e-06,  9.99997715e-01)'

    @test all(m0 .≈ Orient.orient_rot3_icrf_to_mod(j2000s(e)))

    e = Epoch("2025-08-14T11:32:03.804000 TT")
    m0 = DCM(
        9.99980493e-01, -5.72876325e-03, -2.48896620e-03, 
        5.72876341e-03,  9.99983590e-01, -7.06436102e-06, 
        2.48896583e-03, -7.19447528e-06,  9.99996902e-01)'

    @test all(m0 .≈ Orient.orient_rot3_icrf_to_mod(j2000s(e)))

    e = Epoch("2029-07-08T08:17:00.204000 TT")
    m0 = DCM(
        9.99974104e-01, -6.60056979e-03, -2.86769780e-03, 
        6.60056999e-03,  9.99978216e-01, -9.39453283e-06, 
        2.86769734e-03, -9.53415048e-06,  9.99995888e-01)'

    @test all(m0 .≈ Orient.orient_rot3_icrf_to_mod(j2000s(e)))

end