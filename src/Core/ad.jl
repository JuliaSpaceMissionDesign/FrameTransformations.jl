
const TagAD1{T} = ForwardDiff.Tag{Autodiff.JSMDDiffTag, T}
const DualAD1{T} = ForwardDiff.Dual{TagAD1{T}, T, 1}

const TagAD2{T} = ForwardDiff.Tag{Autodiff.JSMDDiffTag, DualAD1{T}}
const DualAD2{T} = ForwardDiff.Dual{TagAD2{T}, DualAD1{T}, 1}

# ------------------------------------------------------------------------------------------
# Points 

const FramePointFunSignature{D, T} = FunctionWrapper{SVector{D, T}, Tuple{T}}
const FramePointFunWrapper{D, T} = FunctionWrappersWrapper{
    Tuple{
        FramePointFunSignature{D, T},
        FramePointFunSignature{D, DualAD1{T}},
        FramePointFunSignature{D, DualAD2{T}}
    }, true
}

function FramePointFunWrapper{D, T}(fun::Function) where {D, T}
    types = (T, DualAD1{T}, DualAD2{T})
    inps = map(x->Tuple{x}, types)
    outs = map(x->SVector{D, x}, types)
    wrps = map(inps, outs) do A, R 
        FunctionWrapper{R, A}(fun)
    end
    return FramePointFunWrapper{D, T}(wrps)
end

# ------------------------------------------------------------------------------------------
# Axes

const FrameAxesFunSignature{O, T} = FunctionWrapper{Rotation{O, T}, Tuple{T}}
const FrameAxesFunWrapper{O, T} = FunctionWrappersWrapper{
    Tuple{
        FrameAxesFunSignature{O, T},
        FrameAxesFunSignature{O, DualAD1{T}},
        FrameAxesFunSignature{O, DualAD2{T}}
    }, true
}

function FrameAxesFunWrapper{O, T}(fun::Function) where {O, T}
    types = (T, DualAD1{T}, DualAD2{T})
    inps = map(x->Tuple{x}, types)
    outs = map(x->Rotation{O, x}, types)
    wrps = map(inps, outs) do A, R 
        FunctionWrapper{R, A}(fun)
    end
    return FrameAxesFunWrapper{O, T}(wrps)
end