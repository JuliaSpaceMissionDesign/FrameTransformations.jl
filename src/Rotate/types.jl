export DCM, EulerAngleAxis, EulerAngles, Quaternion, Rotations

"""
    DCM{T}

Direction Cosine Matrix (DCM) of type `T`, which is a 3x3 static matrix of type
`T`.

### Examples

```julia-repl
julia> DCM(1.0I)
DCM{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0

julia> DCM([1 0 0; 0 -1 0; 0 0 -1])
DCM{Int64}:
 1   0   0
 0  -1   0
 0   0  -1
```
"""
struct DCM{T} <: StaticMatrix{3, 3, T}
    data::NTuple{9, T}

    function DCM(x::NTuple{9, T}) where T
        return new{T}(x)
    end

    function DCM(x::NTuple{9, Any})
        T = StaticArrays.promote_tuple_eltype(x)
        return new{T}(StaticArrays.convert_ntuple(T, x))
    end

    function DCM{T}(x::NTuple{9, T}) where T
        return new{T}(x)
    end

    function DCM{T}(x::NTuple{9, Any}) where T
        return new{T}(StaticArrays.convert_ntuple(T, x))
    end
end

"""
    EulerAngles{T}

The definition of Euler Angles, which is composed of three angles `a1`, `a2`,
and `a3` together with a rotation sequence `seq`.

### Fields

- `a1` -- First rotation (rad).
- `a2` -- Second rotation (rad).
- `a3` -- Third rotation (rad).
- `seq` -- Rotation sequence.

!!! info
    `seq` is provided by a symbol with three characters, each one indicating
    the rotation axis of the corresponding angle, *e.g.* `:ZYX`. The valid
    values for `seq` are:

    - `:XYX`, `:XYZ`, `:XZX`, `:XZY`, `:YXY`, `:YXZ`, `:YZX`, `:YZY`, `:ZXY`,
        `:ZXZ`, `:ZYX`, and `ZYZ`.

### Constructor

    EulerAngles(a1::T1, a2::T2, a3::T3, seq::Symbol = :ZYX) where {T1, T2, T3}

Create a new instance of `EulerAngles` with the angles `a1`, `a2`, and `a3` and
the rotation sequence `seq`.

The type will be inferred from `T1`, `T2`, and `T3`.

If `seq` is not provided, then it defaults to `:ZYX`.

### Examples

```julia-repl
julia> EulerAngles(pi / 2, pi / 4, -pi, :XYZ)
EulerAngles{Float64}:
  R(X) :  1.5707963267948966 rad  ( 90.0°)
  R(Y) :  0.7853981633974483 rad  ( 45.0°)
  R(Z) : -3.141592653589793  rad  (-180.0°)
```
"""
struct EulerAngles{T}
    a1::T
    a2::T
    a3::T
    seq::Symbol
end

function EulerAngles(a1::T1, a2::T2, a3::T3, seq::Symbol = :ZYX) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    return EulerAngles(T(a1), T(a2), T(a3), seq)
end


"""
    EulerAngleAxis{T}

The definition of Euler Angle and Axis to represent a 3D rotation.

### Fields

- `a` -- The Euler angle (rad).
- `v` -- The unitary vector aligned with the Euler axis.

### Constructor

    EulerAngleAxis(a::T, v::AbstractVector{T}) where {T,T}

Create an Euler Angle and Axis representation structure with angle `a` (rad) and
vector `v`.

The vector `v` will not be normalized.

The returned structure type will be selected according to the input types.

### Examples

```julia-repl
julia> EulerAngleAxis(pi / 3, [sqrt(2), sqrt(2), 0])
EulerAngleAxis{Float64}:
  Euler angle:   1.0472 rad ( 60.0000 deg)
   Euler axis: [  1.4142,   1.4142,   0.0000]
```
"""
struct EulerAngleAxis{T}
    a::T
    v::SVector{3, T}
    EulerAngleAxis(a::T, v::SVector{3, T}) where T<:Number = new{T}(a, v)
end

"""
    struct _EulerAngleConversion{R}

This private structure is used only to enable the rotation conversion to Euler
angles using the Julia API.
"""
struct _EulerAngleConversion{R}
end

function EulerAngles(seq::Symbol)
    return _EulerAngleConversion{seq}
end

"""
    Quaternion{T}

The definition of the quaternion.

### Fields

- `q0` -- Quaternion real part.
- `q1` -- X component of the quaternion imaginary part.
- `q2` -- Y component of the quaternion imaginary part.
- `q3` -- Z component of the quaternion imaginary part.

!!! note
    The quaternion `q` in this structure is represented by:

        q = q0 + q1.i + q2.j + q3.k

### Example

```julia-repl
julia> Quaternion(cosd(45), sind(45), 0, 0)
Quaternion{Float64}:
  + 0.7071067811865476 + 0.7071067811865476.i + 0.0.j + 0.0.k
```
"""
struct Quaternion{T}
    q0::T
    q1::T
    q2::T
    q3::T
end

"""
    Rotations

A `Union` of all supported rotation types.
"""
const Rotations = Union{
    DCM,
    EulerAngles,
    EulerAngleAxis,
    Quaternion
}