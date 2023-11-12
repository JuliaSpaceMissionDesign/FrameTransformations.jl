export prepare_eop, read_eop, init_eop, eop_data_filename, IERS_EOP

"""
    EOPData{T}

EOP Data for IAU 2000A.

### Fields
- `filename` : File where the EOP data are stored. 
- `j2000` : Independent valiable - time - in UTC.
- `j2000_TT` : Independent valiable - time - in TT.
- `x, y`: Polar motion with respect to the crust (arcsec).
- `UT1_UTC`: Irregularities of the rotation angle (s).
- `UT1_TT`: Irregularities of the rotation angle (s) w.r.t. TT timescale.
- `LOD`: Length of day offset (ms).
- `dX, dY`: Celestial pole offsets referred to the model IAU2000A (milliarcsec).
"""
mutable struct EOPData{T}
    filename::String
    j2000::Vector{T}
    x::Vector{T}
    y::Vector{T}
    UT1_UTC::Vector{T}
    LOD::Vector{T}
    dX::Vector{T}
    dY::Vector{T}

    # EOP data parametrized by TT epoch
    # These are parsed automatically in eop.jl to allow direct computation of EOP 
    # without performing many transformations from TT/TDB, which are considered to be equal.
    j2000_TT::Vector{T}
    UT1_TT::Vector{T}
end

function EOPData(::Type{T}=Float64) where T 
    return EOPData(
        "", [],[],[],[],[],[],[],
        [],[]
    )
end

function Base.show(io::IO, eop::EOPData{T}) where T 
    if isempty(eop.j2000)
        println(io, "EOPData()")
        return
    end
    println(io, "EOPData(filename=\"$(eop.filename), \"beg=\"$(eop.j2000[1]) (UTC)\", "
        * " end=\"$(eop.j2000[end]) (UTC)\")")
end

"""
    read_iers_eop_finals(filename::AbstractString)

Read IERS EOP C04 files in csv format.


### Returns
- `j2000_utc`: Julian days since J2000 in UTC.
- `j2000_tt`: Julian days since J2000 in TT (Terrestrial Time).
- `x_pole`: Celestial pole offset in the X direction (arcsec).
- `y_pole`: Celestial pole offset in the Y direction (arcsec).
- `ut1_utc`: UT1-UTC time difference (s).
- `ut1_tt`: UT1-TT time difference (s).
- `lod`: Length of Day (s).
- `dX`: Celestial pole offset rate in the X direction (milliarcsec).
- `dY`: Celestial pole offset rate in the Y direction (milliarcsec).

The function reads IERS EOP C04 files in CSV format and extracts relevant Earth Orientation 
Parameters (EOP) data. It then updates predictions, filling missing values with zeros for 
LOD, dX, and dY. 
Finally, the function parametrizes EOP with respect to both UTC and TT time scales for 
convenience.

### References
- https://maia.usno.navy.mil/ser7/readme.finals2000A
- http://hpiers.obspm.fr/eoppc/bul/bulb/explanatory.html
- https://maia.usno.navy.mil
"""
function read_iers_eop_finals(filename::AbstractString)

    data, ~ = readdlm(filename, ';'; header=true)

    last_id_utc = findlast(!isempty, @view(data[:, 15]))
    ut1_utc  = convert(Vector{Float64},  @view(data[1:last_id_utc, 15]))

    mjd = convert(Vector{Float64},  @view(data[1:last_id_utc, 1])) 
    x_pole = convert(Vector{Float64}, @view(data[1:last_id_utc, 6]))
    y_pole = convert(Vector{Float64}, @view(data[1:last_id_utc, 8]))

    lod_raw = @view(data[1:last_id_utc, 17])
    dX_raw = @view(data[1:last_id_utc, 24])
    dY_raw = @view(data[1:last_id_utc, 26])

    # Update predictions
    # Elements that are not present are filled with the last valid element
    last_id_lod = findlast(!isempty, lod_raw) 
    lod_raw[last_id_lod:end] .= lod_raw[last_id_lod-1]
    lod = convert(Vector{Float64}, lod_raw)

    last_id_dXY = findlast(!isempty, dX_raw)
    dX_raw[last_id_dXY:end] .= dX_raw[last_id_dXY-1]
    dY_raw[last_id_dXY:end] .= dY_raw[last_id_dXY-1]
    dX = convert(Vector{Float64}, dX_raw)
    dY = convert(Vector{Float64}, dY_raw)

    # For convenience and to avoid discontinuities, eop are parametrized both in terms of 
    # utc and tt time scales
    j2000_utc = mjd .- Tempo.DMJD
    j2000_tt = map(t -> Tempo.utc2tai(Tempo.DJ2000, t)[2] + Tempo.OFFSET_TAI_TT/Tempo.DAY2SEC, j2000_utc)
    j2000_ut1 = j2000_utc + ut1_utc./Tempo.DAY2SEC  # utc + ut1-utc
    ut1_tt = (j2000_ut1 - j2000_tt) .* Tempo.DAY2SEC 

    return j2000_utc, j2000_tt, x_pole, y_pole, ut1_utc, ut1_tt, lod, dX, dY  
end

"""
    prepare_eop(iers_file::AbstractString, output_filename::AbstractString="iau2000a")  

Prepare Earth Orientation Parameters (EOP) data from IERS EOP C04 files to JSMD's `eop.dat` 
convenience format. The `output_filename` should not include the file extension, which is 
automatically added by this function. 

```@raw julia 
# Save a new file called: test.eop.dat
prepare_eop("input.csv", "test")
```
"""
function prepare_eop(iers_file::AbstractString, output_filename::AbstractString="iau2000a")
    data = hcat(read_iers_eop_finals(iers_file)...)
    writedlm(output_filename * ".eop.dat", data)

    @info "IERS EOP file '$(iers_file)' converted to '$(output_filename).eop.dat'."

    nothing
end

"""
    read_eop(filename)  

Read Earth Orientation Parameters (EOP) from JSMD `.eop.dat` file.  
"""
function read_eop(filename)
    !endswith(filename, "eop.dat") && throw(
        ArgumentError("EOP reader support only '.eop.dat' files! Please prepare " * 
                       "the data with \'prepare_eop\' and retry."))

    # Load the file
    data = readdlm(filename; header=false)
    
    j2000_utc = @view(data[:, 1])
    j2000_tt = @view(data[:, 2])
    x_pole = @view(data[:, 3])
    y_pole = @view(data[:, 4])
    ut1_utc = @view(data[:, 5])
    ut1_tt = @view(data[:, 6])
    lod = @view(data[:, 7])
    dX = @view(data[:, 8])
    dY = @view(data[:, 9])

    return j2000_utc, j2000_tt, x_pole, y_pole, ut1_utc, ut1_tt, lod, dX, dY  
end

"""
    set_eop_data(filename)  

Set Earth Orientation Parameters (EOP) to be used for frames transformations from JSMD 
`.eop.dat` file.  
"""
function set_eop_data(filename::AbstractString)

    oldfile = IERS_EOP_DATA.filename
    j2000_utc, j2000_tt, x_pole, y_pole, ut1_utc, ut1_tt, lod, dX, dY = read_eop(filename)

    if (!isempty(IERS_EOP_DATA.j2000))
        @warn "Existing EOP data from \'$(oldfile)\' will be overwritten by \'$(filename)\'."
    end

    IERS_EOP_DATA.filename = filename
    IERS_EOP_DATA.j2000 = j2000_utc
    IERS_EOP_DATA.x = x_pole
    IERS_EOP_DATA.y = y_pole
    IERS_EOP_DATA.UT1_UTC = ut1_utc
    IERS_EOP_DATA.LOD = lod
    IERS_EOP_DATA.dX = dX
    IERS_EOP_DATA.dY = dY
    IERS_EOP_DATA.j2000_TT = j2000_tt
    IERS_EOP_DATA.UT1_TT = ut1_tt

    nothing
end

"""
    eop_data_filename()  

Get loaded Earth Orientation Parameters (EOP) filename.
"""
function eop_data_filename()
    if isempty(IERS_EOP_DATA.filename) 
        throw(ErrorException("Unable to retrieve filename, no EOP data has been loaded."))
    end
    return IERS_EOP_DATA.filename
end

"""
    EOPInterpolator{T}

EOP Data for IAU 2000A.

### Fields
- `init`: A flag indicating whether the EOPInterpolator has been initialized.
- `j2000`: Independent variable (time), in UTC.
- `j2000_TT`: Independent variable (time), in TT.
- `x, y`: Polar motion with respect to the crust (arcsec).
- `UT1_UTC`: Irregularities of the rotation angle (s).
- `UT1_TT`: Irregularities of the rotation angle (s) with respect to TT timescale.
- `LOD`: Length of day offset (ms).
- `dX, dY`: Celestial pole offsets referred to the model IAU2000A (milliarcsec).
- `x_TT, y_TT`: Polar motion with respect to the crust (arcsec) parametrized by TT epoch.
- `UT1_TT`: Irregularities of the rotation angle (s) parametrized by TT epoch.
- `LOD_TT`: Length of day offset (ms) parametrized by TT epoch.
- `dX_TT, dY_TT`: Celestial pole offsets referred to the model IAU2000A (milliarcsec) parametrized by TT epoch.
"""
mutable struct EOPInterpolator{T<:AbstractInterpolationMethod}
    init::Bool

    x::T
    y::T
    UT1_UTC::T
    LOD::T
    dX::T
    dY::T

    # EOP data parametrized by TT epoch
    # These are parsed automatically in eop.jl to allow direct computation of EOP 
    # without performing many transformations from TT/TDB, which are considered to be equal.
    x_TT::T
    y_TT::T
    UT1_TT::T
    LOD_TT::T
    dX_TT::T
    dY_TT::T
end 

function Base.show(io::IO, s::EOPInterpolator{T}) where T 
    println(io, "EOPInterpolator(init=$(s.init))")
end

function _eop_spline_initializer(::Type{T}=Float64) where T
    return InterpAkima([-100005.0, -100004.0, -100003.0, -100002.0, -100001.0, -100000.0], zeros(6))
end

function EOPInterpolator(::Type{T}=Float64) where T
    EOPInterpolator(
        false, 
        _eop_spline_initializer(), _eop_spline_initializer(), _eop_spline_initializer(),
        _eop_spline_initializer(), _eop_spline_initializer(), _eop_spline_initializer(),
        _eop_spline_initializer(), _eop_spline_initializer(), _eop_spline_initializer(),
        _eop_spline_initializer(), _eop_spline_initializer(), _eop_spline_initializer()
    )
end

"""
    init_eop(filename)  

Initialize Earth Orientation Parameters (EOP) from file.  

!!! warn
    This function must be called to initialize the EOP data used by frames, in case 
    Earth-associated frames are used.  

!!! warn 
    This function accept only `.eop.dat` files. Please use [`prepare_eop`](@ref) to transform 
    IERS EOP files in this format.
"""
function init_eop(
    filename::AbstractString, ::Type{INTERP} = InterpAkima) where {INTERP <: AbstractInterpolationMethod}
    
    # Set eop data to gather
    set_eop_data(filename)

    # Initialize and set interpolators
    IERS_EOP.init = true
    IERS_EOP.x = INTERP(IERS_EOP_DATA.j2000, IERS_EOP_DATA.x)
    IERS_EOP.y = INTERP(IERS_EOP_DATA.j2000, IERS_EOP_DATA.y)
    IERS_EOP.UT1_UTC = INTERP(IERS_EOP_DATA.j2000, IERS_EOP_DATA.UT1_UTC)
    IERS_EOP.LOD = INTERP(IERS_EOP_DATA.j2000, IERS_EOP_DATA.LOD)
    IERS_EOP.dX = INTERP(IERS_EOP_DATA.j2000, IERS_EOP_DATA.dX)
    IERS_EOP.dY = INTERP(IERS_EOP_DATA.j2000, IERS_EOP_DATA.dY)
    
    IERS_EOP.x_TT = INTERP(IERS_EOP_DATA.j2000_TT, IERS_EOP_DATA.x)
    IERS_EOP.y_TT = INTERP(IERS_EOP_DATA.j2000_TT, IERS_EOP_DATA.y)
    IERS_EOP.UT1_TT = INTERP(IERS_EOP_DATA.j2000_TT, IERS_EOP_DATA.UT1_TT)
    IERS_EOP.LOD_TT = INTERP(IERS_EOP_DATA.j2000_TT, IERS_EOP_DATA.LOD)
    IERS_EOP.dX_TT = INTERP(IERS_EOP_DATA.j2000_TT, IERS_EOP_DATA.dX)
    IERS_EOP.dY_TT = INTERP(IERS_EOP_DATA.j2000_TT, IERS_EOP_DATA.dY)

    @info "EOP initialized from file '$(filename)'."

    nothing 

end

"""
    IERS_EOP_DATA

Earth Orientation Parameters Data: x/y pole, UT1-UTC, LOD, dX, dY (smoothed values at 1-day 
intervals) with respect to IAU 2006/2000A precession-nutation model and consistent with ITRF2014.

See also: [`EOPData`](@ref)
"""
const IERS_EOP_DATA = EOPData();

"""
    IERS_EOP 

Earth Orientation Parameters interpolators: x/y pole, UT1-UTC, LOD, dX, dY (smoothed values 
at 1-day intervals) with respect to IAU 2006/2000A precession-nutation model and consistent 
with ITRF2014.

See also: [`EOPInterpolator`](@ref)
"""
const IERS_EOP = EOPInterpolator();

"""
    offset_utc2ut1(seconds)

Return the offset between `UTC` and `UT1` in seconds.
"""
@inline function offset_utc2ut1(seconds)
    utc = seconds / 86400.0
    !IERS_EOP.init && throw(
        ErrorException(
            "EOP not initialized. Please run 'init_eop' before using this function."
        ))
    return interpolate(IERS_EOP.UT1_UTC, utc)
end