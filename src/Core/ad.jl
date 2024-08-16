const TagAD1{T} = ForwardDiff.Tag{JSMDDiffTag, T}
const DualAD1{T} = ForwardDiff.Dual{TagAD1{T}, T, 1}

# ------------------------------------------------------------------------------------------
# Points 

const FramePointFunSignature{O, T} = FunctionWrapper{Translation{O, T}, Tuple{T}}
const FramePointFunWrapper{O, T} = FunctionWrappersWrapper{
    Tuple{
        FramePointFunSignature{O, T},
        FramePointFunSignature{O, DualAD1{T}}
    }, true
}

function FramePointFunWrapper{O, T}(fun::Function) where {O, T}
    types = (T, DualAD1{T})
    inps = map(x->Tuple{x}, types)
    outs = map(x->Translation{O, x}, types)
    wrps = map(inps, outs) do A, R 
        FunctionWrapper{R, A}(fun)
    end
    return FramePointFunWrapper{O, T}(wrps)
end


# ------------------------------------------------------------------------------------------
# Axes

const FrameAxesFunSignature{O, T} = FunctionWrapper{Rotation{O, T}, Tuple{T}}
const FrameAxesFunWrapper{O, T} = FunctionWrappersWrapper{
    Tuple{
        FrameAxesFunSignature{O, T},
        FrameAxesFunSignature{O, DualAD1{T}}
    }, true
}

function FrameAxesFunWrapper{O, T}(fun::Function) where {O, T}
    types = (T, DualAD1{T})
    inps = map(x->Tuple{x}, types)
    outs = map(x->Rotation{O, x}, types)
    wrps = map(inps, outs) do A, R 
        FunctionWrapper{R, A}(fun)
    end
    return FrameAxesFunWrapper{O, T}(wrps)
end