struct Date 
    year::Int
    month::Int 
    day::Int
    function Date(y::N, m::N, d::N) where N
        return new(y, m, d)
    end
end