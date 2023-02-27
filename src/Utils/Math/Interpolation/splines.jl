export InterpCubicSplines

struct InterpCubicSplines{N, S, T} <: AbstractInterpolationMethod
    n::Int
    x::Vector{N}
    y::Matrix{N}
    c::Array{N, S}
    s::Int
end

@inbounds function InterpCubicSplines(x::Vector{N}, y::Matrix{N}, 
    type::Symbol=:Natural) where N 
    n = length(x)
    n < 4 && throw(ArgumentError("At least four points are needed."))
    row, col = size(y)
    col != n && throw(ArgumentError("`x` and `y` must have the same length."))
    c = Vector{N}()
    for i in 1:row
        append!(c, _create_splines(Val(type), x, y[i, :]))
    end
    return InterpCubicSplines{N, row, type}(
        n, collect(x), collect(y), N.(reshape(c, (4, col-1, row))), row)
end

@inbounds function CubicSplines(x::Vector{N}, y::Matrix{N}, 
    type::Symbol=:Natural) where N 
    n = length(x)
    n < 4 && throw(ArgumentError("At least four points are needed."))
    row, col = size(y)
    col != n && throw(ArgumentError("`x` and `y` must have the same length."))
    c = Vector{N}()
    for i in 1:row
        append!(c, _create_splines(Val(type), x, y[i, :]))
    end
    return CubicSplines{N, row, type}(
        n, collect(x), collect(y), N.(reshape(c, (4, col-1, row))), row)
end

@inbounds function _create_splines(::Val{:NotAKnot}, x::Vector{N}, y::Vector{N}) where N
    n = length(x)
    dl = similar(y, n-1)
    d = similar(y, n)
    du = similar(y, n-1)
    b = similar(y, n)
    c = similar(y, 4, n-1)
    @views begin
        dx = diff(x)
        slope = diff(y) ./ dx

        d[2:end-1] .= 2 .* (dx[1:end-1] .+ dx[2:end])
        du[2:end] .= dx[1:end-1]
        dl[1:end-1] .= dx[2:end]
        b[2:end-1] .= 3 * (dx[2:end] .* slope[1:end-1]
            .+ dx[1:end-1] .* slope[2:end])

        # Not-a-knot boundary conditions
        d[1] = dx[2]
        du[1] = x[3] - x[1]
        δ = x[3] - x[1]
        b[1] = ((dx[1] + 2δ) * dx[2] * slope[1]
            + dx[1]^2 * slope[2]) / δ
        d[end] = dx[end-1]
        dl[end] = x[end] - x[end-2]
        δ = x[end] - x[end-2]
        b[end] = ((dx[end]^2 * slope[end-1]
                + (2δ + dx[end]) * dx[end-1] * slope[end]) / δ)

        A = Tridiagonal(dl, d, du)
        s = A \ b
        t = (s[1:end-1] .+ s[2:end] .- 2 .* slope) ./ dx

        c[1, :] .= y[1:end-1]
        c[2, :] .= s[1:end-1]
        c[3, :] .= (slope .- s[1:end-1]) ./ dx .- t
        c[4, :] .= t ./ dx
    end
    return c 
end

@inbounds function _create_splines(::Val{:Natural}, x::Vector{N}, y::Vector{N}) where N
    n = length(x)
    dl = similar(y, n-1)
    d = similar(y, n)
    du = similar(y, n-1)
    b = similar(y, n)
    c = similar(y, 4, n-1)
    @views begin
        dx = diff(x)
        slope = diff(y) ./ dx

        d[2:end-1] .= 2 .* (dx[1:end-1] .+ dx[2:end])
        du[2:end] .= dx[1:end-1]
        dl[1:end-1] .= dx[2:end]
        b[2:end-1] .= 3 * (dx[2:end] .* slope[1:end-1]
            .+ dx[1:end-1] .* slope[2:end])

        # Natural boundary conditions
        d[1] = 1
        b[1] = 0
        d[end] = 1
        b[end] = 0

        A = Tridiagonal(dl, d, du)
        s = A \ b
        t = (s[1:end-1] .+ s[2:end] .- 2 .* slope) ./ dx

        c[1, :] .= y[1:end-1]
        c[2, :] .= s[1:end-1]
        c[3, :] .= (slope .- s[1:end-1]) ./ dx .- t
        c[4, :] .= t ./ dx
    end
    return c 
end

@inbounds function _create_splines(::Val{:Periodic}, x::Vector{N}, y::Vector{N}) where N
    n = length(x)
    dl = similar(y, n-1)
    d = similar(y, n)
    du = similar(y, n-1)
    b = similar(y, n)
    c = similar(y, 4, n-1)
    @views begin
        dx = diff(x)
        slope = diff(y) ./ dx

        # Periodic boundary conditions
        slope[1] = (y[2] - y[end]) / (x[2] - x[end])
        slope[end] = (y[1] - y[end-1]) / (x[1] - x[end-1])

        d[1:end-1] .= 2 .* (dx[1:end-1] .+ dx[2:end])
        d[end] = 2*(dx[end-1] + dx[1])
        du[2:end] .= dx[1:end-1]
        du[1] = dx[end]
        dl[1:end-1] .= dx[2:end]
        dl[end] = dx[1]
        b[1:end-1] .= 3 * (dx[2:end] .* slope[1:end-1]
            .+ dx[1:end-1] .* slope[2:end])
        b[end] = 3*(dx[1]*slope[end]+dx[end]*slope[1])

        A = Tridiagonal(dl, d, du)
        s = A \ b
        t = (s[1:end-1] .+ s[2:end] .- 2 .* slope) ./ dx

        c[1, :] .= y
        c[2, :] .= s
        c[3, :] .= (slope .- s[1:end-1]) ./ dx .- t
        c[4, :] .= t ./ dx
    end
    return c 
end

@inbounds function _create_splines(::Val{:Quadratic}, x::Vector{N}, y::Vector{N}) where N
    n = length(x)
    dl = similar(y, n-1)
    d = similar(y, n)
    du = similar(y, n-1)
    b = similar(y, n)
    c = similar(y, 4, n-1)

    @views begin
        dx = diff(x)
        slope = diff(y) ./ dx

        d[2:end-1] .= 2 .* (dx[1:end-1] .+ dx[2:end])
        du[2:end] .= dx[1:end-1]
        dl[1:end-1] .= dx[2:end]
        b[2:end-1] .= 3 * (dx[2:end] .* slope[1:end-1]
            .+ dx[1:end-1] .* slope[2:end])

        # Quadratic boundary conditions
        d[1] = 2*dx[1] + dx[2]
        du[1] = dx[1]
        b[1] = 3*dx[1]*slope[1] + 3*dx[2]*slope[1] - 3*dx[1]*slope[2]
        
        d[end] = 2*dx[end-1] + dx[end]
        dl[end] = dx[end-1]
        b[end] = 3*dx[end-1]*slope[end-1] + 3*dx[end]*slope[end-1] - 3*dx[end]*slope[end-2]

        A = Tridiagonal(dl, d, du)
        s = A \ b
        t = (s[1:end-1] .+ s[2:end] .- 2 .* slope) ./ dx

        c[1, :] .= y[1:end-1]
        c[2, :] .= s[1:end-1]
        c[3, :] .= (slope .- s[1:end-1]) ./ dx .- t
        c[4, :] .= t ./ dx
    end

    return c 
end

function interpolate(cs::InterpCubicSplines{N, S, T}, x) where {N, S, T}
    j = searchsortedfirst(cs.x, x)
    if j == 1
        y0 = @view cs.y[1:end, 1]
        return SVector{S}(yi for yi in y0)
    elseif j == cs.n
        yn = @view cs.y[1:end, end]
        return SVector{S}(yi for yi in yn)
    end
    dx = x - cs.x[j-1]

    return SVector{S}(
        cs.c[1, j-1, i] + dx * (cs.c[2, j-1, i] + dx * (cs.c[3, j-1, i] + dx * cs.c[4, j-1, i])) for i in 1:S
    )
end