function orient_body_angles(b, T::N) where {N<:AbstractFloat}
    orient_right_ascension(b, T) + π/2, π/2 - orient_declination(b, T), mod2pi(orient_rotation_angle(b, T))
end

function orient_body_rates(b, T::N) where {N<:AbstractFloat}
    orient_right_ascension_rate(b, T), -orient_declination_rate(b, T), orient_rotation_rate(b, T)
end