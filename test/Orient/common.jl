
function _test_matrix(e, g, fun)
    v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

    ft = fun(j2000s(e))
    err = g*ft' - I(3)
    @test norm(err) ≤ 1e-8

    for i in 1:3 
        v = zeros(3)
        v[i] = 1.0
        @test v2as(g * v, ft * v) ≤ 1e-8
        v[i] = -1.0
        @test v2as(g * v, ft * v) ≤ 1e-8
    end
end

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
