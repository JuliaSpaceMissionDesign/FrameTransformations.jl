
# ==========================================================================================
# Compiled fast-path: opt-in API that extracts raw closures from FunctionWrappers,
# giving users zero-overhead, AD-transparent callables for the hot path.
# ==========================================================================================

"""
    _extract_raw_fn(fww::FunctionWrappersWrapper)

Extract the original Julia closure from a `FunctionWrappersWrapper`. The raw closure
preserves its concrete type, enabling the compiler to inline and differentiate through it —
unlike `FunctionWrapper` which erases the type behind a C function pointer.

!!! note "Internal layout dependency"
    This accesses `fww.fw[1].obj[]` — the first `FunctionWrapper` in the dispatch tuple
    (the `Float64` signature), then dereferences the `Base.RefValue` holding the closure.
    This depends on the internal layout of FunctionWrappers.jl (v1.x) and
    FunctionWrappersWrappers.jl (v0.1.x). The compat bounds in Project.toml must be kept
    tight to guard against silent breakage if these packages change internals.
"""
_extract_raw_fn(fww::FunctionWrappersWrapper) = fww.fw[1].obj[]

# ------------------------------------------------------------------------------------------
# CompiledRotation
# ------------------------------------------------------------------------------------------

"""
    CompiledRotation{O, Inv, F}

A zero-overhead, AD-transparent callable that computes a `Rotation{O}` between two axes.
Created via [`compile_rotation`](@ref). The type parameter `F` captures the concrete closure
type, allowing the compiler to inline the call. The `Inv` parameter encodes whether the
rotation should be inverted, eliminating the runtime branch entirely.

### Usage
```julia
cr = compile_rotation(fr, :ICRF, :BODY)
R  = cr(t)   # ::Rotation{O} — no FunctionWrapper overhead
```
"""
struct CompiledRotation{O,Inv,F}
    fun::F
end

function CompiledRotation{O,Inv}(fun::F) where {O,Inv,F}
    return CompiledRotation{O,Inv,F}(fun)
end

@inline function (cr::CompiledRotation{O,false})(t::Number) where {O}
    return Rotation{O}(cr.fun(t))
end

@inline function (cr::CompiledRotation{O,true})(t::Number) where {O}
    return inv(Rotation{O}(cr.fun(t)))
end

# ------------------------------------------------------------------------------------------
# CompiledTranslation
# ------------------------------------------------------------------------------------------

"""
    CompiledTranslation{O, F}

A zero-overhead, AD-transparent callable that computes a `Translation{O}` (as `SVector`)
between two points in a specified axes frame. Created via [`compile_translation`](@ref).

### Usage
```julia
ct = compile_translation(fr, :Earth, :Moon, :ICRF)
v  = ct(t)   # ::SVector — no FunctionWrapper overhead
```
"""
struct CompiledTranslation{O,F}
    fun::F
end

function CompiledTranslation{O}(fun::F) where {O,F}
    return CompiledTranslation{O,F}(fun)
end

@inline function (ct::CompiledTranslation{O})(t::Number) where {O}
    return ct.fun(t)
end

# ------------------------------------------------------------------------------------------
# CompiledDirection
# ------------------------------------------------------------------------------------------

"""
    CompiledDirection{O, F}

A zero-overhead, AD-transparent callable that computes a direction vector at a given time.
Created via [`compile_direction`](@ref).

### Usage
```julia
cd = compile_direction(fr, :SunDir, :ICRF)
d  = cd(t)   # ::SVector — no FunctionWrapper overhead
```
"""
struct CompiledDirection{O,F}
    fun::F
end

function CompiledDirection{O}(fun::F) where {O,F}
    return CompiledDirection{O,F}(fun)
end

@inline function (cd::CompiledDirection{O})(t::Number) where {O}
    return cd.fun(t)
end

# ------------------------------------------------------------------------------------------
# compile_rotation
# ------------------------------------------------------------------------------------------

"""
    compile_rotation(fr::FrameSystem, from, to)
    compile_rotation(fr::FrameSystem, from, to, ::Val{N})

Compile a zero-overhead rotation callable between axes `from` and `to`. The returned
[`CompiledRotation`](@ref) bypasses `FunctionWrapper` dispatch entirely by extracting
the raw closure and preserving its concrete type.

By default the compiled callable uses the frame system's maximum order `O`. Pass
`Val{N}()` where `N ≤ O` to extract only the `N`-th order closure, avoiding unnecessary
derivative computation in the hot path.

Supports both direct parent-child pairs and multi-hop paths. For multi-hop paths,
the construction loop is type-unstable (closures are composed via `let`), but the
resulting callable is fully typed — so `cr(t)` in a hot loop is zero-overhead.

!!! warning
    The compiled callable captures a snapshot of the current frame graph topology.
    If axes are added after compilation, the callable becomes stale and must be
    recompiled.
"""
function compile_rotation(fr::FrameSystem{O,T}, from, to) where {O,T}
    return compile_rotation(fr, from, to, Val(O))
end

function compile_rotation(fr::FrameSystem{O,T}, from, to, ::Val{N}) where {O,T,N}
    N > O && throw(
        ArgumentError("requested order $N exceeds frame system order $O.")
    )

    fromid = axes_id(fr, from)
    toid = axes_id(fr, to)

    fromid == toid && return CompiledRotation{N,false}(_ -> one(T) * I)

    nodes = _get_axes_nodes(fr, fromid, toid)
    isnothing(nodes) && throw(
        ErrorException("no path between axes $fromid and $toid in the frame system.")
    )

    return _compile_rotation(Val(N), nodes)
end

function _compile_rotation_pair(::Val{N}, from::FrameAxesNode, to::FrameAxesNode) where {N}
    if from.id == to.parentid
        raw_fn = _extract_raw_fn(to.f[N])
        return CompiledRotation{N,false}(raw_fn)
    else
        raw_fn = _extract_raw_fn(from.f[N])
        return CompiledRotation{N,true}(raw_fn)
    end
end

function _compile_rotation(::Val{N}, nodes::Vector{<:FrameAxesNode}) where {N}
    cr = _compile_rotation_pair(Val(N), nodes[1], nodes[2])

    for i in 3:length(nodes)
        cr_next = _compile_rotation_pair(Val(N), nodes[i-1], nodes[i])
        cr = let inner = cr, outer = cr_next
            CompiledRotation{N,false}(t -> outer(t) * inner(t))
        end
    end

    return cr
end

# ------------------------------------------------------------------------------------------
# compile_translation
# ------------------------------------------------------------------------------------------

"""
    compile_translation(fr::FrameSystem, from, to, axes)
    compile_translation(fr::FrameSystem, from, to, axes, ::Val{N})

Compile a zero-overhead translation callable between points `from` and `to`, expressed in
the given `axes` frame. The returned [`CompiledTranslation`](@ref) bypasses `FunctionWrapper`
dispatch for both the point functions and any axis rotations needed along the path.

By default the compiled callable uses the frame system's maximum order `O`. Pass
`Val{N}()` where `N ≤ O` to extract only the `N`-th order closure.

!!! warning
    The compiled callable captures a snapshot of the current frame graph topology.
    If points or axes are added after compilation, the callable becomes stale and
    must be recompiled.
"""
function compile_translation(fr::FrameSystem{O,T}, from, to, axes) where {O,T}
    return compile_translation(fr, from, to, axes, Val(O))
end

function compile_translation(fr::FrameSystem{O,T}, from, to, axes, ::Val{N}) where {O,T,N}
    N > O && throw(
        ArgumentError("requested order $N exceeds frame system order $O.")
    )

    fromid = point_id(fr, from)
    toid = point_id(fr, to)
    axid = axes_id(fr, axes)

    if fromid == toid
        return CompiledTranslation{N}(_ -> @SVector zeros(T, 3 * N))
    end

    nodes = _get_points_nodes(fr, fromid, toid)
    isnothing(nodes) && throw(
        ErrorException("no path between points $fromid and $toid in the frame system.")
    )

    return _compile_translation(Val(N), fr, nodes, axid)
end

function _compile_point_pair(::Val{N}, from::FramePointNode, to::FramePointNode) where {N}
    if from.id == to.parentid
        return to.axesid, _extract_raw_fn(to.f[N]), false
    else
        return from.axesid, _extract_raw_fn(from.f[N]), true
    end
end

function _compile_translation(::Val{N}, fr::FrameSystem, nodes::Vector{<:FramePointNode}, axes::Int) where {N}
    if length(nodes) == 2
        axid, raw_fn, inv_flag = _compile_point_pair(Val(N), nodes[1], nodes[2])
        if axid != axes
            cr = compile_rotation(fr, axid, axes, Val(N))
            if inv_flag
                return CompiledTranslation{N}(let _fn = raw_fn, _cr = cr
                    t -> SVector(_cr(t) * (-Translation{N}(_fn(t))))
                end)
            else
                return CompiledTranslation{N}(let _fn = raw_fn, _cr = cr
                    t -> SVector(_cr(t) * Translation{N}(_fn(t)))
                end)
            end
        else
            if inv_flag
                return CompiledTranslation{N}(let _fn = raw_fn
                    t -> SVector(-Translation{N}(_fn(t)))
                end)
            else
                return CompiledTranslation{N}(let _fn = raw_fn
                    t -> SVector(Translation{N}(_fn(t)))
                end)
            end
        end
    end

    return _compile_translation_forward(Val(N), fr, nodes, axes)
end

function _compile_translation_forward(::Val{N}, fr::FrameSystem, nodes::Vector{<:FramePointNode}, out_axes::Int) where {N}
    axid1, raw1, inv1 = _compile_point_pair(Val(N), nodes[1], nodes[2])

    compiled_fun = let _fn = raw1, _inv = inv1
        if _inv
            t -> (axid1, -Translation{N}(_fn(t)))
        else
            t -> (axid1, Translation{N}(_fn(t)))
        end
    end

    prev_axid = axid1

    for i in 3:length(nodes)
        axid_i, raw_i, inv_i = _compile_point_pair(Val(N), nodes[i-1], nodes[i])

        if axid_i != prev_axid
            cr = compile_rotation(fr, prev_axid, axid_i, Val(N))
            compiled_fun = let _prev = compiled_fun, _fn = raw_i, _inv = inv_i, _cr = cr, _axid = axid_i
                function (t)
                    _, tr = _prev(t)
                    tr_rotated = _cr(t) * tr
                    tr_hop = _inv ? -Translation{N}(_fn(t)) : Translation{N}(_fn(t))
                    return (_axid, tr_rotated + tr_hop)
                end
            end
        else
            compiled_fun = let _prev = compiled_fun, _fn = raw_i, _inv = inv_i, _axid = axid_i
                function (t)
                    _, tr = _prev(t)
                    tr_hop = _inv ? -Translation{N}(_fn(t)) : Translation{N}(_fn(t))
                    return (_axid, tr + tr_hop)
                end
            end
        end

        prev_axid = axid_i
    end

    if prev_axid != out_axes
        cr_final = compile_rotation(fr, prev_axid, out_axes, Val(N))
        return CompiledTranslation{N}(let _inner = compiled_fun, _cr = cr_final
            function (t)
                _, tr = _inner(t)
                return SVector(_cr(t) * tr)
            end
        end)
    else
        return CompiledTranslation{N}(let _inner = compiled_fun
            function (t)
                _, tr = _inner(t)
                return SVector(tr)
            end
        end)
    end
end

# ------------------------------------------------------------------------------------------
# compile_direction
# ------------------------------------------------------------------------------------------

"""
    compile_direction(fr::FrameSystem, name::Symbol, axes)
    compile_direction(fr::FrameSystem, name::Symbol, axes, ::Val{N})

Compile a zero-overhead direction callable for the direction `name`, expressed in the
given `axes` frame. The returned [`CompiledDirection`](@ref) bypasses `FunctionWrapper`
dispatch.

By default the compiled callable uses the frame system's maximum order `O`. Pass
`Val{N}()` where `N ≤ O` to extract only the `N`-th order closure.

!!! warning
    The compiled callable captures a snapshot of the current frame graph topology.
    If axes are added after compilation, the callable becomes stale and must be
    recompiled.
"""
function compile_direction(fr::FrameSystem{O}, name::Symbol, axes) where {O}
    return compile_direction(fr, name, axes, Val(O))
end

function compile_direction(fr::FrameSystem{O}, name::Symbol, axes, ::Val{N}) where {O,N}
    N > O && throw(
        ArgumentError("requested order $N exceeds frame system order $O.")
    )

    if !has_direction(fr, name)
        throw(
            ErrorException("No direction with name $name registered in the frame system.")
        )
    end

    node = directions(fr)[name]
    raw_fn = _extract_raw_fn(node.f[N])
    thisaxid = node.axesid
    axid = axes_id(fr, axes)

    if thisaxid != axid
        cr = compile_rotation(fr, thisaxid, axid, Val(N))
        return CompiledDirection{N}(let _fn = raw_fn, _cr = cr
            function (t)
                stv = Translation{N}(_fn(t))
                return SVector(_cr(t) * stv)
            end
        end)
    else
        return CompiledDirection{N}(let _fn = raw_fn
            function (t)
                return SVector(Translation{N}(_fn(t)))
            end
        end)
    end
end
