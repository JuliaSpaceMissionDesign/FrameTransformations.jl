@testset "IERS Transformations" verbose=true begin 

    @testset "Fundamental Arguments" begin 

        # Testing IERS 2010 conventions Fundamental Arguments 
        # Test values are in radians and taken from ERFA 
        # Absolute tolerance of 1e-13 guarantees an error below 0.02 μas

        # Testing at t = 1. to check correct coefficients!

        t = [0.06567, 1.]
        fa = [Orient.FundamentalArguments(t[1]), Orient.FundamentalArguments(t[2])]

        atol, rtol = 1e-13, 1e-12 

        # Testing Delaunay's Arguments 
        @testset "Delaunay arguments 2003" begin

            # Testing Mean Anomaly of the Moon 
            @test fa[1].Mₐ ≈ 2.663600612437975    atol=atol rtol=rtol
            @test fa[2].Mₐ ≈ 5.826604253498457    atol=atol rtol=rtol

            # Testing Mean Anomaly of the Sun 
            @test fa[1].Sₐ ≈ 3.518352361195792    atol=atol rtol=rtol
            @test fa[2].Sₐ ≈ 6.223481898965822    atol=atol rtol=rtol

            # Testing Mean Argument of Latitude of the Moon 
            @test fa[1].uₘ ≈ 2.533320307830868    atol=atol rtol=rtol
            @test fa[2].uₘ ≈ 3.059317938336433    atol=atol rtol=rtol

            # Testing mean elongation of the moon from the sun 
            @test fa[1].Dₛ ≈ 3.236084178870326e-1 atol=atol rtol=rtol
            @test fa[2].Dₛ ≈ 4.275356347495070    atol=atol rtol=rtol

            # Testing Mean Longitude of the Moon
            @test fa[1].Ωₘ ≈ -3.438585492123376e-2 + 2π atol=atol rtol=rtol
            @test fa[2].Ωₘ ≈ -1.586439578169723e-1 + 2π atol=atol rtol=rtol

        end

        @testset "Planetary Arguments" begin

            # Testing Mean Longitude of Mercury
            @test fa[1].λ_Me ≈ 6.075865478867676 atol=atol rtol=rtol
            @test fa[2].λ_Me ≈ 5.671020519871966 atol=atol rtol=rtol

            # Testing Mean Longitude of Venus
            @test fa[1].λ_Ve ≈ 1.131754499992205    atol=atol rtol=rtol
            @test fa[2].λ_Ve ≈ 3.454962478273771e-1 atol=atol rtol=rtol

            # Testing Mean Longitude of Earth
            @test fa[1].λ_Ea ≈ 5.315317577813381 atol=atol rtol=rtol
            @test fa[2].λ_Ea ≈ 1.742524595141305 atol=atol rtol=rtol

            # Testing Mean Longitude of Mars
            @test fa[1].λ_Ma ≈ 3.008541490420559 atol=atol rtol=rtol
            @test fa[2].λ_Ma ≈ 9.727169953023775e-1 atol=atol rtol=rtol

            # Testing Mean Longitude of Jupiter
            @test fa[1].λ_Ju ≈ 4.078027048663447 atol=atol rtol=rtol
            @test fa[2].λ_Ju ≈ 3.303160303663311 atol=atol rtol=rtol

            # Testing Mean Longitude of Saturn
            @test fa[1].λ_Sa ≈ 2.274751979272320 atol=atol rtol=rtol
            @test fa[2].λ_Sa ≈ 3.354371331461241 atol=atol rtol=rtol

            # Testing Mean Longitude of Uranus
            @test fa[1].λ_Ur ≈ 5.972384629789489 atol=atol rtol=rtol
            @test fa[2].λ_Ur ≈ 3.930831143408273e-1 atol=atol rtol=rtol

            # Testing Mean Longitude of Neptune
            @test fa[1].λ_Ne ≈ 5.562305932034747 atol=atol rtol=rtol
            @test fa[2].λ_Ne ≈ 2.842004543620414 atol=atol rtol=rtol

            # Testing General Accumulated Precession 
            @test fa[1].pₐ ≈ 1.601172753812795e-3 atol=atol rtol=rtol
            @test fa[2].pₐ ≈ 2.438713691000000e-2 atol=atol rtol=rtol

        end 

        # Testing the truncated expressions of the Delunary Arguments 
        fab = [Orient.FundamentalArguments(t[1], iau2006b), 
            Orient.FundamentalArguments(t[2], iau2006b)];

        @testset "Delaunay arguments IAU2000B" begin 

            # Testing Mean Anomaly of the Moon 
            @test fab[1].Mₐ ≈ 2.663599945842266      atol=atol rtol=rtol
            @test fab[2].Mₐ ≈ 5.826449449627781      atol=atol rtol=rtol

            # Testing Mean Anomaly of the Sun 
            @test fab[1].Sₐ ≈ 3.518352372771510          atol=atol rtol=rtol
            @test fab[2].Sₐ ≈ 6.223484580361229          atol=atol rtol=rtol

            # Testing Mean Argument of Latitude of the Moon 
            @test fab[1].uₘ ≈ 2.533320574432192          atol=atol rtol=rtol
            @test fab[2].uₘ ≈ 3.059379762905646          atol=atol rtol=rtol

            # Testing mean elongation of the moon from the sun 
            @test fab[1].Dₛ ≈ 3.236085510636535e-1       atol=atol rtol=rtol
            @test fab[2].Dₛ ≈ 4.275387201216140          atol=atol rtol=rtol

            # Testing Mean Longitude of the Moon
            @test fab[1].Ωₘ ≈ -3.438601115926876e-2 + 2π  atol=atol rtol=rtol 
            @test fab[2].Ωₘ ≈ -1.586802211172697e-1 + 2π  atol=atol rtol=rtol 

            for f in [:λ_Me, :λ_Ve, :λ_Ea, :λ_Ma, :λ_Ju, :λ_Sa, :λ_Ur, :λ_Ne, :pₐ]
                for i in eachindex(fab)
                    @test getproperty(fab[i], f) == getproperty(fa[i], f)
                end
            end

        end 

    end;


    @testset "Nutation" begin 

        jdₜ = 2400000.5 + 53750.892855; 
        ERFA_DJ00 = 2451545.
        ERFA_DJC = 36525.

        t = ((jdₜ - ERFA_DJ00)) / ERFA_DJC;

        r2a = 180/π*3600 

        # Current Nutation series are taken from ERFA 
        @testset "IAU 2006 A/B Nutation Models" begin

            # Testing IAU 2000A Nutation model from ERFA 
            # must be precise up to 1 μas
            atol, rtol = 1e-7, 1e-7

            fa = Orient.FundamentalArguments(t);
            Δψ, Δϵ = Orient.nutation00(iau2006a, t, fa) .* r2a

            @test Δψ ≈ -1.071332651430803 atol=atol rtol=rtol
            @test Δϵ ≈  8.656842463521295 atol=atol rtol=rtol

            Δψ, Δϵ = Orient.orient_nutation(iau2006a, t) .* r2a
            @test Δψ ≈ -1.071332974891340 atol=atol rtol=rtol
            @test Δϵ ≈  8.656841011106838 atol=atol rtol=rtol
            
            # Testing IAU 2000B Nutation model from ERFA 
            # Must be precise up to n 1 mas 
            atol, rtol = 1e-12, 1e-12
            Δψ, Δϵ = Orient.orient_nutation(iau2006b, t) .* r2a

            @test Δψ ≈ -1.071752875965778 atol=atol rtol=rtol
            @test Δϵ ≈  8.656781912467901 atol=atol rtol=rtol
        end  
    end

    @testset "Precession" begin 

        jdₜ = 2400000.5 + 53750.892855; 
        ERFA_DJ00 = 2451545.
        ERFA_DJC = 36525.

        t = ((jdₜ - ERFA_DJ00)) / ERFA_DJC;

        r2a = 180/π*3600 
        atol, rtol = 1e-12, 1e-12

        # Testing Fukushima-Williams angles in μas 
        # Values taken from ERFA, error must be below 1 μas to be acceptable.

        fw = [Orient.fw_angles(iau2006a, t).* r2a, Orient.fw_angles(iau2006a, 0).* r2a] 

        @testset "Fukushima-Williams angles" begin 
            # Testing γ
            @test fw[1][1] ≈ 5.865586624386230e-1 atol=atol rtol=rtol
            @test fw[2][1] ≈ -5.292800000000158e-2 atol=atol rtol=rtol

            # Testing ϕ 
            @test fw[1][2] ≈ 8.437858525780872e4 atol=atol rtol=rtol
            @test fw[2][2] ≈ 8.438141281900251e4 atol=atol rtol=rtol

            # Testing longitude ψ
            @test fw[1][3] ≈ 3.043272121523480e2 atol=atol rtol=rtol
            @test fw[2][3] ≈ -4.177500000000124e-2 atol=atol rtol=rtol

            # Testing obliquity ϵ
            @test fw[1][4] ≈ 8.43785766962175e4 atol=atol rtol=rtol
            @test fw[2][4] ≈ 8.438140600000251e4 atol=atol rtol=rtol
        end

        @testset "Rotation Matrices" begin 

            γ, ϕ, ψ, ϵ = Orient.fw_angles(iau2006a, t); 
            FW = Orient.fw_matrix(γ, ϕ, ψ, ϵ)

            # Check the matrix originating from the FW Precession Angles 
            @test FW[1, 1] ≈ 9.999989154136046e-1  atol=atol rtol=rtol
            @test FW[1, 2] ≈ -1.350835287774888e-3 atol=atol rtol=rtol
            @test FW[1, 3] ≈ -5.868693548984236e-4 atol=atol rtol=rtol
            @test FW[2, 1] ≈  1.350835311644812e-3 atol=atol rtol=rtol
            @test FW[2, 2] ≈  9.999990876215008e-1 atol=atol rtol=rtol
            @test FW[2, 3] ≈ -3.557088188443913e-7 atol=atol rtol=rtol
            @test FW[3, 1] ≈  5.868692999554672e-4 atol=atol rtol=rtol
            @test FW[3, 2] ≈ -4.370554148591665e-7 atol=atol rtol=rtol
            @test FW[3, 3] ≈  9.999998277921019e-1 atol=atol rtol=rtol

            # Check the whole Assembly process 
            FW2 = Orient.orient_precession_bias(iau2006a, t);

            @test FW2[1, 1] ≈ 9.999989154136046e-1  atol=atol rtol=rtol
            @test FW2[1, 2] ≈ -1.350835287774888e-3 atol=atol rtol=rtol
            @test FW2[1, 3] ≈ -5.868693548984236e-4 atol=atol rtol=rtol
            @test FW2[2, 1] ≈  1.350835311644812e-3 atol=atol rtol=rtol
            @test FW2[2, 2] ≈  9.999990876215008e-1 atol=atol rtol=rtol
            @test FW2[2, 3] ≈ -3.557088188443913e-7 atol=atol rtol=rtol
            @test FW2[3, 1] ≈  5.868692999554672e-4 atol=atol rtol=rtol
            @test FW2[3, 2] ≈ -4.370554148591665e-7 atol=atol rtol=rtol
            @test FW2[3, 3] ≈  9.999998277921019e-1 atol=atol rtol=rtol

        end

    end

    @testset "ITRF to GCRF Routines" begin

        r2a = 180/π*3600 

        jdₜ = 2400000.5 + 53750.892855; 
        ERFA_DJ00 = 2451545.
        ERFA_DJC = 36525.

        t = ((jdₜ - ERFA_DJ00)) / ERFA_DJC;

        atol, rtol = 1e-12, 1e-12

        @testset "Polar Motion" begin 
            # Testing TIO Locator Position [in arcseconds] 
            sp = Orient.tio_locator(t) .* r2a
            @test sp ≈ -2.839163974949028e-6 atol=atol rtol=rtol

            # Polar Motion Matrix 
            xₚ, yₚ = 1.857, 0.123;  
            W = Orient.polar_motion(xₚ, yₚ, t)

            @test W[1, 1] ≈ -2.823123663590987e-1   atol=atol rtol=rtol
            @test W[1, 2] ≈  1.176993683001583e-1   atol=atol rtol=rtol
            @test W[1, 3] ≈ -9.520748849236963e-1   atol=atol rtol=rtol
            @test W[2, 1] ≈  3.885932432356596e-12  atol=atol rtol=rtol
            @test W[2, 2] ≈  9.924450321335735e-1   atol=atol rtol=rtol
            @test W[2, 3] ≈  1.226900900374203e-1   atol=atol rtol=rtol
            @test W[3, 1] ≈  9.593225358557600e-1   atol=atol rtol=rtol
            @test W[3, 2] ≈  3.463692964357531e-2   atol=atol rtol=rtol
            @test W[3, 3] ≈ -2.801795055034182e-1   atol=atol rtol=rtol
        end

        @testset "Earth Rotation Angle" begin 
            
            # Testing Earth Rotation Angle [radians]
            # Note in this case the input time should be expressed in UT1
            era = Orient.earth_rotation_angle(jdₜ); 
            @test era ≈ 1.335810933027503 atol=atol rtol=rtol
        end

        @testset "Precession-Nutation" begin 

            # Testing CIP Coordinates from FW angles 
            # Error should be below 1μas
            x, y = Orient.fw2xy(0.1, 0.7, 1.1, -2.5);

            @test x ≈ -5.614951267014441e-1 atol=atol rtol=rtol        
            @test y ≈  2.536947155044642e-1 atol=atol rtol=rtol

            atol, rtol = 1e-7, 1e-7
            # Testing X, Y coordinates including IAU2000A nutation series 
            # Error must be below 1 μas
            x, y = Orient.cip_coords(iau2006a, t) .* r2a;
            @test x ≈ 1.206359974143122e2 atol=atol rtol=rtol
            @test y ≈ 8.567258642517157   atol=atol rtol=rtol

            # Testing CIO Locator [in radians]
            # Error should be below 1 μas 
            s = Orient.cio_locator(iau2006a, t, 0.125, 0.745); 
            @test s ≈ -4.656250032319267e-2 atol=1e-12 rtol=1e-12
            
            atol, rtol = 1e-12, 1e-12

            # Testing IAU 2006A CIP Motion Rotation Matrix 
            Q = Orient.cip_motion(iau2006a, t); 
            @test Q[1, 1] ≈  9.999998289694806e-1   atol=atol rtol=rtol
            @test Q[1, 2] ≈ -2.461548571919270e-8   atol=atol rtol=rtol
            @test Q[1, 3] ≈  5.848598198075143e-4   atol=atol rtol=rtol
            @test Q[2, 1] ≈  3.231916262391721e-10  atol=atol rtol=rtol
            @test Q[2, 2] ≈  9.999999991374118e-1   atol=atol rtol=rtol
            @test Q[2, 3] ≈  4.153524199496106e-5   atol=atol rtol=rtol
            @test Q[3, 1] ≈ -5.848598203254312e-4   atol=atol rtol=rtol
            @test Q[3, 2] ≈ -4.153523470214526e-5   atol=atol rtol=rtol
            @test Q[3, 3] ≈  9.999998281068927e-1   atol=atol rtol=rtol
            
            # Testing IAU2000B CIP Motion Rotation Matrix 
            # Tolerances are 1000x the previous because this model is 1000x less 
            # precise than the IAU2006A 
            atol, rtol = 1e-9, 1e-9
            Q = Orient.cip_motion(iau2006b, t); 
            @test Q[1, 1] ≈  9.999998289699551e-1   atol=atol rtol=rtol
            @test Q[1, 2] ≈ -2.461541308285131e-8   atol=atol rtol=rtol
            @test Q[1, 3] ≈  5.848590084916441e-4   atol=atol rtol=rtol
            @test Q[2, 1] ≈  3.232319550905416e-10  atol=atol rtol=rtol
            @test Q[2, 2] ≈  9.999999991374174e-1   atol=atol rtol=rtol
            @test Q[2, 3] ≈  4.153510643734409e-5   atol=atol rtol=rtol
            @test Q[3, 1] ≈ -5.848590090095588e-4   atol=atol rtol=rtol
            @test Q[3, 2] ≈ -4.153509914454784e-5   atol=atol rtol=rtol
            @test Q[3, 3] ≈  9.999998281073728e-1   atol=atol rtol=rtol

        end

    end

end;