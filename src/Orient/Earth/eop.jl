export get_iers_eop, IERS_EOP

"""
    get_iers_eop(; force_download = false)

Download and parse the IERS EOP C04 data. 

The files are downloaded using the `RemoteFile` package with weekly updates. Hence, if one 
desires to force a download before the scheduled time, then set the keyword `force_download`  
to `true`.

!!! note
    The files will be downloaded from the default URL. If the user want to use
    another one, then use the specialized function [`get_iers_eop_IAU2000A`](@ref)

    See also: [`get_iers_eop_IAU2000A`](@ref)

### Returns
A structure [`EOPData`](@ref) with the interpolations of the EOP parameters. Notice that the
interpolation indexing is set to Julian days since J2000.

"""
function get_iers_eop(; force_download = false)
    return get_iers_eop_IAU2000A(force_download = force_download)
end

"""
    get_iers_eop_iau_2000A(url::String="https://datacenter.iers.org/data/csv/finals2000A.all.csv"; 
        force_download = false)

Get the IERS EOP C04 IAU2000A data from the URL `url`.

If `url` is omitted, then it defaults to https://datacenter.iers.org/data/csv/finals2000A.all.csv.
The file is downloaded using the `RemoteFile` package with weekly updates. Hence, if one desires 
to force a download before the scheduled time, then set the keyword `force_download` to `true`.

!!! note
    The interpolation of every field in [`EOPData`](@ref) between two
    points in the grid is linear. If extrapolation is needed, then if will use
    the nearest value (flat extrapolation).

See also: [`get_iers_eop`](@ref)

### Returns
The structure `EOPData` with the interpolations of the EOP parameters. Notice that the 
interpolation indexing is set to Julian days since J2000.
"""
function get_iers_eop_IAU2000A(
    url::String = "https://datacenter.iers.org/data/csv/finals2000A.all.csv";
    force_download = false
)
    @RemoteFile(
        _eop_iau2000A,
        url,
        file="eop_iau2000a.txt",
        dir=joinpath(@__DIR__, "..", "..", "..", "ext"),
        updates=:fridays,
        failed=:warn
    )

    # Download the data
    download(_eop_iau2000A; force = force_download, force_update = true)

    # Parse the data removing the header.
    eop, ~ = readdlm(path(_eop_iau2000A), ';'; header = true)

    # Obtain the last available index of the field.
    last_id = findlast(!isempty, eop[:, 11])
    last_id === nothing && (last_id = length(eop[:, 11]))
    ut1_utc::Vector{Float64} = Vector{Float64}(eop[1:last_id, 11])

    # Create the EOP Data structure by creating the interpolations:
    # - The interpolation will be linear between two points in the grid.
    # - The extrapolation will be flat, considering the nearest point.
    j2000_utc = Vector{Float64}(eop[1:last_id, 1] .+ 2400000.5 .- DJ2000)
    j2000_tt = [Tempo.utc2tai(DJ2000, utci)[2] for utci in j2000_utc] 
        .- Tempo.OFFSET_TAI_TT ./ Tempo.DAY2SEC

    j2000_ut1 = j2000_utc + Vector{Float64}(ut1_utc) ./ Tempo.DAY2SEC  # utc + ut1-utc
    ut1_tt = (j2000_ut1 - j2000_tt) .* Tempo.DAY2SEC .- Tempo.OFFSET_TAI_TT

    return EOPData(
        _create_iers_eop_interpolation(j2000_utc, eop[:, 6]),
        _create_iers_eop_interpolation(j2000_utc, eop[:, 8]),
        _create_iers_eop_interpolation(j2000_utc, eop[:, 11]),
        _create_iers_eop_interpolation(j2000_utc, eop[:, 13]), # This value is in ms 
        _create_iers_eop_interpolation(j2000_utc, eop[:, 20]), # This value is in mas 
        _create_iers_eop_interpolation(j2000_utc, eop[:, 22]), # This value is in mas

        _create_iers_eop_interpolation(j2000_tt, eop[:, 6]),
        _create_iers_eop_interpolation(j2000_tt, eop[:, 8]),
        _create_iers_eop_interpolation(j2000_tt, ut1_tt),
        _create_iers_eop_interpolation(j2000_tt, eop[:, 13]), # This value is in ms 
        _create_iers_eop_interpolation(j2000_tt, eop[:, 20]), # This value is in mas
        _create_iers_eop_interpolation(j2000_tt, eop[:, 22]), # This value is in mas
    )

end

function Base.show(io::IO, eop::EOPData{T}) where T
    # Check if IO has support for colors.
    println(io, " ")
    println(io, "  EOPData ", "│ ", "Timespan (UTC)")
    println(io, " ─────────┼──────────────────────────────────────────────────")
    println(io, "        x ", "│ ", _get_iers_eop_timespan(eop.x))
    println(io, "        y ", "│ ", _get_iers_eop_timespan(eop.y))
    println(io, "  UT1-UTC ", "│ ", _get_iers_eop_timespan(eop.UT1_UTC))
    println(io, "      LOD ", "│ ", _get_iers_eop_timespan(eop.LOD))
    println(io, "       dX ", "│ ", _get_iers_eop_timespan(eop.dX))
    println(io, "       dY ", "│ ", _get_iers_eop_timespan(eop.dY))

    return nothing
end

# Get timespan
function _get_iers_eop_timespan(itp::InterpAkima)
    str = string(DateTime(first(itp.x)*86400)) * " - " *
          string(DateTime(last(itp.x)*86400))
    return str
end

# Create the interpolation object for the `knots` and `field` from IERS.
function _create_iers_eop_interpolation(
    knots::AbstractVector,
    field::AbstractVector
)
    # Obtain the last available index of the field.
    last_id = findlast(!isempty, field)
    last_id === nothing && (last_id = length(field))

    # Convert the field to a `Vector{Float64}`.
    field_float::Vector{Float64} = Vector{Float64}(field[1:last_id])

    # Create the interpolation object.
    return InterpAkima(knots[1:last_id], field_float)

end

"""
    IERS_EOP 

Earth orientation parameters: x/y pole, UT1-UTC, LOD, dX, dY (smoothed values at 1-day intervals) 
with respect to IAU 2006/2000A precession-nutation model and consistent with ITRF2014.

EOP 14 C04 is updated two times per week.

Here the files are downloaded using the `RemoteFile` package with weekly updates.

See also: [`get_iers_eop`](@ref)
"""
const IERS_EOP = get_iers_eop()



# Adds UT1 Timescale to Tempo!
# It is defined here because referring the variable IERS_EOP inside Tempo would lead to 
# allocations! 

"""
    offset_utc2ut1(seconds)

Return the offset between [`UTC`](@ref) and [`UT1`](@ref) in seconds.
"""
@inline function offset_utc2ut1(seconds)
    utc = seconds/86400.0
    return interpolate(IERS_EOP.UT1_UTC, utc)
end

Tempo.add_timescale(Tempo.TIMESCALES, Tempo.UT1, offset_utc2ut1, parent=Tempo.UTC)