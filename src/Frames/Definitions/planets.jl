using Basic.Bodies: NAIFId

abstract type BodyCentricInertialTrueOfDate <: AbstractBodyCentricInertialFrame end

abstract type BodyCentricInertial2000 <: AbstractBodyCentricInertialFrame end

function BodyCentricInertialTrueOfDate(id::NAIFId, name::Symbol) 
    valId = id.valId
    tname = Symbol("BodyCentricInertial", "_", name)
    sname = Symbol("BCI", "_", uppercase(String(name)))
    return """
        # ------------------------------------------------------------------------------
        # Frame definition: $tname
        # ------------------------------------------------------------------------------

        # Type
        struct $tname <: BodyCentricInertialTrueOfDate
            vid::Val{Int}
            function $tname()
                new($valId)
            end
        end

        # Singleton
        $sname = $tname()

        # Transformations
        function Frames.Rotation(
            origin::InternationalCelestialReferenceFrame,
            target::$(tname), 
            e::Epoch
        )
            t = j2000c(convert(TT, e))
            α = π/2 + orient_right_ascension(vid, t)  
            δ = π/2 - orient_declination(vid, t)
            Rmat = angle_to_dcm(α, δ, :ZX)

            Frames.Rotation(origin, target, Rmat)
        end
        function Frames.Rotation(
            origin::$(tname),
            target::InternationalCelestialReferenceFrame,
            e::Epoch
        )
            inv(Frames.Rotation(target, origin, e))
        end

        # Connect to frame graph
        connect!(ICRF, $sname)

        # Export type and singleton
        export $tname, $sname
    """
end

function BodyCentricInertial2000(id::NAIFId, name::Symbol)
    valId = id.valId
    tname = Symbol("BodyCentricInertial2000", "_", name)
    sname = Symbol("BCI2000", "_", uppercase(String(name)))
    rname = Symbol("ICRF2",uppercase(String(name)),"_", "BCI2000")
    return """
        # ------------------------------------------------------------------------------
        # Frame definition: $tname
        # ------------------------------------------------------------------------------

        # Type
        struct $tname <: BodyCentricInertial2000
            vid::Val{Int}
            function $tname()
                new($valId)
            end
        end

        # Singleton
        $sname = $tname()

        # Declare function to build the BCI2000
        function __build_bci2000rot(::$(valId))
            α = π/2 + orient_right_ascension(vid, 0.0)  
            δ = π/2 - orient_declination(vid, 0.0)
            angle_to_dcm(α, δ, :ZX)
        end
        const $rname = __build_bci2000rot(valId)

        # Transformations
        function Frames.Rotation(
            origin::InternationalCelestialReferenceFrame,
            target::$(tname), 
            e::Epoch
        )   
            Frames.Rotation(origin, target, $rname)
        end
        function Frames.Rotation(
            origin::$(tname),
            target::InternationalCelestialReferenceFrame,
            e::Epoch
        )
            inv(Frames.Rotation(target, origin, e))
        end

        # Connect to frame graph
        connect!(ICRF, $sname)

        # Export type and singleton
        export $tname, $sname
    """
end