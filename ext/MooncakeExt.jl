module MooncakeExt

using FrameTransformations: FrameSystem,
    vector3, vector6, vector9, vector12,
    rotation3, rotation6, rotation9, rotation12,
    direction3, direction6, direction9, direction12

using FunctionWrappers: FunctionWrapper
using FunctionWrappersWrappers: FunctionWrappersWrapper

import Mooncake
using Mooncake: @from_chainrules, DefaultCtx

# ==========================================================================================
# Declare FrameSystem and its internal opaque types as non-differentiable
# ==========================================================================================

# FunctionWrapper and FunctionWrappersWrapper are compiled C function pointers — opaque
# to reverse-mode AD. Declaring them as NoTangent prevents Mooncake from trying to
# construct tangent types for them (which would fail).
Mooncake.tangent_type(::Type{<:FunctionWrapper}) = Mooncake.NoTangent
Mooncake.tangent_type(::Type{<:FunctionWrappersWrapper}) = Mooncake.NoTangent

# FrameSystem is a structural container (graph of axes/points with function wrappers).
# It is never differentiated through — only `t` is differentiable.
Mooncake.tangent_type(::Type{<:FrameSystem}) = Mooncake.NoTangent

# ==========================================================================================
# Register ChainRulesCore rrules with Mooncake
# ==========================================================================================

# Vector functions: vector3, vector6, vector9, vector12
@from_chainrules DefaultCtx Tuple{typeof(vector3), FrameSystem, Any, Any, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(vector6), FrameSystem, Any, Any, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(vector9), FrameSystem, Any, Any, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(vector12), FrameSystem, Any, Any, Any, T} where {T<:Number}

# Rotation functions: rotation3, rotation6, rotation9, rotation12
@from_chainrules DefaultCtx Tuple{typeof(rotation3), FrameSystem, Any, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(rotation6), FrameSystem, Any, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(rotation9), FrameSystem, Any, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(rotation12), FrameSystem, Any, Any, T} where {T<:Number}

# Direction functions: direction3, direction6, direction9, direction12
@from_chainrules DefaultCtx Tuple{typeof(direction3), FrameSystem, Symbol, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(direction6), FrameSystem, Symbol, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(direction9), FrameSystem, Symbol, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(direction12), FrameSystem, Symbol, Any, T} where {T<:Number}

end # module
