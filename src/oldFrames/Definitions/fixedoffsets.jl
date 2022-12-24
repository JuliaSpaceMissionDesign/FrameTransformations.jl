export @frame_fixed_rotation_euler, @frame_fixed_rotation_quaternion

using Random: randstring

"""
    frame_fixed_rotation_euler(graph, parent, axis, angles)

Create a fixed offset rotation reference frame.

### Inputs 
- `graph` -- graph where to connect the frame 
- `parent` -- parent frame type -- shall have an overload with no arguments
- `axis` -- axis or axes of rotation 
- `angles` -- angles of rotation
"""
macro frame_fixed_rotation_euler(graph, parent, axis, angles)
    s = randstring(20)
    R = Symbol("ROTMAT", "_", "FIXROTEULER", s)
    uniquename = Symbol(:FixedRotationEulerFrame, s)
    return quote 
        # type 
        struct $uniquename <: FixedRotationFrame end
        # rotation matrix
        const $R = angle_to_dcm($angles..., $axis)
        # rotations
        function Rotation(origin::$parent, target::$uniquename, e::Epoch)
            new(origin, target, $R)
        end
        function Rotation(origin::$uniquename, target::$parent, e::Epoch)
            new(origin, target, transpose($(R)))
        end
        # connect to graph
        connect!($graph, $(parent)(), $(uniquename)())
    end
end

"""
    frame_fixed_rotation_quaternion(graph, parent, axis, angles)

Create a fixed offset rotation reference frame.

### Inputs 
- `graph` -- graph where to connect the frame 
- `parent` -- parent frame type -- shall have an overload with no arguments
- `quat` -- quaternion
"""
macro frame_fixed_rotation_quaternion(graph, parent, quat)
    s = randstring(20)
    R = Symbol("ROTMAT", "_", "FIXROTQUAT", s)
    uniquename = Symbol(:FixedRotationQuaterionFrame, s)
    return quote 
        # type 
        struct $uniquename <: FixedRotationFrame end
        # rotation matrix
        const $R = quat_to_angle($angles..., $axis)
        # rotations
        function Rotation(origin::$parent, target::$uniquename, e::Epoch)
            new(origin, target, $R)
        end
        function Rotation(origin::$uniquename, target::$parent, e::Epoch)
            new(origin, target, transpose($(R)))
        end
        # connect to graph
        connect!($graph, $(parent)(), $(uniquename)())
    end
end