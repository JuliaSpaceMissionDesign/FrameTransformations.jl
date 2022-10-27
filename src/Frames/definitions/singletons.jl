const FRAMES2SINGL = (
    (:ICRF, :InternationalCelestialReferenceFrame),
    (:MEOD, :MeanEclipticOfDate),
    (:MEJ2000, :MeanEclipticJ2000),
)

for (s, n) in FRAMES2SINGL
    @eval begin
        """
            $s

        Singleton instance of the [`$n`](@ref).
        """
        const $s = $n()
        export $s, $n
    end
end