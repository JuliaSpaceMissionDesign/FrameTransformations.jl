

struct TestC{T} 
    f::FunctionWrapper{T, Tuple{T, Int}}
end

function build_test(fn)

    dfun(t) = cos(fn(t))
    TestC{Float64}((t, x) -> dfun(t))

end

x = build_test(sin)
xr = x.f(π/6, 1)

y = build_test(cos)
yr = y.f(π/6, 1)