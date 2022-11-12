export @frame_bci_tod_iau

macro frame_bci_tod_iau(body)
    tname = Symbol(body, :BodyCentricInertialFrame)
    sname = Symbol("BCI_TOD", "_", uppercase(String(body)))

    return quote 
        # Frame type
        struct $tname  <: BodyCentricInertialTrueOfDateFrame end

        function Frames.build(::Type{$tname}, args...)
            return $tname()
        end

        const $sname = $tname()  # Singleton

        # Transformation
        function Rotation(
            origin::InternationalCelestialReferenceFrame, target::$tname,
            e::Epoch
        )
            t = j2000c(convert(TT, e))
            α = π/2 + orient_right_ascension($body, t)  
            δ = π/2 - orient_declination($body, t)
            R = angle_to_dcm(-α, -δ, :ZX)
            Rotation(origin, target, R)
        end
        function Rotation(
            origin::$tname, target::InternationalCelestialReferenceFrame, 
            e::Epoch
        )
            inv(Rotation(target, origin, e))
        end

        connect!(FRAMES, $sname) # Connect to frame graph
        export $tname, $sname # Export type and singleton
    end
end


macro frame_bci_2000_iau(body)
    tname = Symbol(body, :BodyCentricInertial2000Frame)
    sname = Symbol("BCI_2000", "_", uppercase(String(body)))
    rname = Symbol("_ROTMAT_", "ICRF2",uppercase(String(name)),"_", "BCI2000")
    α = π/2 + orient_right_ascension($body, 0.0)  
    δ = π/2 - orient_declination($body, 0.0)
    R = angle_to_dcm(-α, -δ, :ZX)

    return quote
        struct $tname <: BodyCentricInertial2000Frame end 

        function Frames.build(::Type{$tname}, args...)
            return $tname()
        end

        const $sname = $tname()  # Singleton
        const $rname = $R

        # Transformation
        function Rotation(
            origin::InternationalCelestialReferenceFrame, target::$tname,
            e::Epoch
        )
            Rotation(origin, target, R)
        end
        function Rotation(
            origin::$tname, target::InternationalCelestialReferenceFrame, 
            e::Epoch
        )
            inv(Rotation(target, origin, e))
        end

        connect!(FRAMES, $sname) # Connect to frame graph
        export $tname, $sname # Export type and singleton
    end
end