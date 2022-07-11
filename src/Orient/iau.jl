function iau_angles(body::AbstractBody, ep::Float64)
    α = body_right_ascention(body, ep)
    δ = body_declination(body, ep)
    θ = body_rotation_angle(body, ep)
    return EulerAngles(α, δ, θ)
end