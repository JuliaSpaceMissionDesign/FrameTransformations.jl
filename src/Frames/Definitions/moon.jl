export add_axes_pa440!, 
       add_axes_pa421!, 
       add_axes_me421! 


"""
	add_axes_pa440!(frames, axes, parent) 

Add `axes` as a set of ephemeris axes representing the DE440 Moon's Principal Axes (PA) to 
`frames`. The libration angles are extracted from the ephemeris kernels loaded within `frames`, 
an error is thrown if such orientation data is not available. 

The `parent` set of axes must be named `ICRF` or have ID = 1 (i.e., the International 
Celestial Reference Frame), otherwise and error is thrown.

!!! warning
	To properly read the ephemeris kernels, the ID associated to the input `axes` must match 
	NAIF's FRAME ID for the Moon PA DE440 axes (31008).

### See also 
See also [`Orient.AXESID_MOONPA_DE440`](@ref), [`Orient.orient_rot3_icrf_to_pa440`](@ref), 
[`add_axes_pa421!`](@ref), and [`add_axes_me421!`](@ref)
"""
function add_axes_pa440!(frames::FrameSystem, axes::AbstractFrameAxes, 
			parent::AbstractFrameAxes)

	if axes_name(parent) != :ICRF && axes_id(parent) != Orient.AXESID_ICRF
		throw(
			ArgumentError("The DE440 Moon Principal Axes (PA) can only be defined "*
				"w.r.t. the International Celestial Reference Frame (ICRF)."
			)
		)
	end

	axesid = axes_id(axes)
	if axesid != Orient.AXESID_MOONPA_DE440 
		throw(
			ArgumentError("$(axes_name(axes)) is associated to an ID that is not the standard 
				Moon PA DE440 ID. Found $axesid, expected ($(Orient.AXESID_MOONPA_DE440))."
			)
		)
	end

	add_axes_ephemeris!(frames, axes, :ZXZ)

end


"""
	add_axes_pa421!(frames, axes, parent) 

Add `axes` as a set of ephemeris axes representing the DE421 Moon's Principal Axes (PA) to 
`frames`. The libration angles are extracted from the ephemeris kernels loaded within `frames`, 
an error is thrown if such orientation data is not available. 

The `parent` set of axes must be named `ICRF` or have ID = 1 (i.e., the International 
Celestial Reference Frame), otherwise and error is thrown.

!!! warning
	To properly read the ephemeris kernels, the ID associated to the input `axes` must match 
	NAIF's FRAME ID for the Moon PA DE421 axes (31006).

### See also 
	See also [`Orient.AXESID_MOONPA_DE421`](@ref), [`Orient.orient_icrf_to_pa421`](@ref), 
	[`add_axes_pa440!`](@ref), and [`add_axes_me421!`](@ref)
"""
function add_axes_pa421!(frames::FrameSystem, axes::AbstractFrameAxes, 
			parent::AbstractFrameAxes)

	if axes_name(parent) != :ICRF && axes_id(parent) != Orient.AXESID_ICRF
		throw(
			ArgumentError("The DE421 Moon Principal Axes (PA) can only be defined "*
				"w.r.t. the International Celestial Reference Frame (ICRF)."
			)
		)
	end

	axesid = axes_id(axes)
	if axesid != Orient.AXESID_MOONPA_DE421
		throw(
			ArgumentError("$(axes_name(axes)) is associated to an ID that is not the standard 
				Moon PA DE421 ID. Found $axesid, expected ($(Orient.AXESID_MOONPA_DE421))."
			)
		)
	end

	add_axes_ephemeris!(frames, axes, :ZXZ)

end


"""
	add_axes_mer421!(frames, axes, parent) 

Add `axes` as fixed offset axes representing the DE421 Moon's Mean Earth/Mean Rotation (MER) 
to `frames`.

The `parent` set of axes can be either the DE440 Principal Axes (PA440) or the DE421 
Principal Axes (PA421), otherwise an error is thrown. Depending on that, the relative axes 
orientation will be automatically selected by this function. 

### See also 
	See also [`add_axes_pa440!`](@ref), and [`add_axes_pa421!`](@ref), 
	[`Orient.DCM_MOONPA421_TO_MER421`](@ref) and [`Orient.DCM_MOONPA421_TO_MER421`](@ref), 
	
"""
function add_axes_me421!(frames::FrameSystem, axes::AbstractFrameAxes, 
			parent::AbstractFrameAxes)

	pid = axes_id(parent) 

	if pid == Orient.AXESID_MOONPA_DE421
		dcm = Orient.DCM_MOONPA421_TO_MER421

	elseif pid == Orient.AXESID_MOONPA_DE440
		dcm = Orient.DCM_MOONPA440_TO_MER421
		
	else 
		pname = axes_name(parent)
		throw(ArgumentError(
			"The DE421 Mean Earth/Mean Rotation (MER) axes cannot be defined w.r.t. "*
			"$pname with ID $pid. Only the DE440 or DE421 Moon Principal Axes are "*
			"accepted as parent axes."
		))
	end

	add_axes_fixedoffset!(frames, axes, parent, dcm)
end