
export skew 


skew(a) = SMatrix{3, 3}(0, a[3], -a[2], -a[3], 0, a[1], a[2], -a[1], 0)


function angle_to_δ³dcm(θ, rot_seq::Symbol)

    s, c = sincos(θ[1])

    a, b = (s, c).*θ[2].*θ[3].*3
    d, e = (s, c).*θ[2]^3
    f, g = (s, c).*θ[4]

    if rot_seq == :X 
        return DCM(0, 0, 0, 0, -f-b+d, -g+a+e, 0, g-a-e, -f-b+d) 

    elseif rot_seq == :Y 
        return DCM(-f-b+d, 0, g-a-e, 0, 0, 0, -g+a+e, 0, -f-b+d) 

    elseif rot_seq == :Z 
        return DCM(-f-b+d, -g+a+e, 0, g-a-e, -f-b+d, 0, 0, 0, 0) 

    else 
        throw(ArgumentError("rot_seq must be :X, :Y or :Z"))
    end

end

function angle_to_δ²dcm(θ, rot_seq::Symbol)

    s, c = sincos(θ[1])

    a, b = (s, c).*θ[2]^2
    d, e = (s, c).*θ[3]

    if rot_seq == :X 
        return DCM(0, 0, 0, 0, -d-b, a-e, 0, e-a, -d-b )

    elseif rot_seq == :Y 
        return DCM(-d-b, 0, e-a, 0, 0, 0, a-e, 0, -d-b)

    elseif rot_seq == :Z
        return DCM(-d-b, a-e, 0, e-a, -d-b,  0, 0, 0, 0) 

    else
        throw(ArgumentError("rot_seq must be :X, :Y, or :Z"))
    end

end


function angle_to_δdcm(θ, rot_seq::Symbol)
    
    s, c = sincos(θ[1]).*θ[2]

    if rot_seq == :X 
        return DCM(0, 0, 0, 0, -s, -c, 0, c, -s)

    elseif rot_seq == :Y 
        return DCM(-s, 0, c, 0, 0, 0, -c, 0, -s)

    elseif rot_seq == :Z
        return DCM(-s, -c, 0, c, -s, 0, 0, 0, 0)

    else
        throw(ArgumentError("rot_seq must be :X, :Y, or :Z"))
    end

end

using ReferenceFrameRotations



function angle_to_δdcm(θ, ϕ, rot_seq::Symbol)

    s, c = sincos(ϕ[1])
    b, a = sincos(θ[1])

    δϕ, δθ = θ[2], ϕ[2]

    if rot_seq == :YX
        return DCM(
            -δθ*b, δϕ*c*b + δθ*s*a, δθ*c*a - δϕ*s*b, 
            0, -δϕ*s, -δϕ*c, 
            -δθ*a, δϕ*c*a - δθ*s*b, -δϕ*s*a - δθ*c*b
        )
    
    elseif rot_seq == :ZX
        return DCM(
            -δθ*b, δϕ*s*b - δθ*c*a, δϕ*c*b + δθ*s*a, 
            δθ*a, δθ*c*b - δϕ*s*a, -δϕ*c*a - δθ*s*b, 
            0, δϕ*c, -δϕ*s
        )

    elseif rot_seq == :YX 

    elseif rot_seq == :YZ 

    elseif rot_seq == :ZX 

    elseif rot_seq == :ZY 
    
    else 
        throw(ArgumentError(
            "The rotation sequence :$rot_seq is not valid."
        ))
    end

end

function angle_to_δdcm(θ, ϕ, rot_seq::Symbol)

    s, c = sincos(θ[1])
    b, a = sincos(ϕ[1])

    δϕ, δθ = ϕ[2], θ[2]

    if rot_seq == :XY
        return DCM(
            -b*δϕ, 0, a*δϕ, 
            a*s*δϕ + b*c*δθ, -s*δθ, b*s*δϕ-a*c*δθ, 
            b*s*δθ -a*c*δϕ, c*δθ, -b*c*δϕ-a*s*δθ
        )
        
    elseif rot_seq == :XZ
        return DCM(
            -b*δϕ, -a*δϕ, 0, 
            a*c*δϕ -b*s*δθ, -b*c*δϕ-a*s*δθ, -c*δθ, 
            a*s*δϕ+b*c*δθ, a*c*δθ-b*s*δϕ, -s*δθ
        )

    elseif rot_seq == :YX 
        return DCM(
            -δθ*s, δθ*c*b + δϕ*s*a, δθ*c*a - δϕ*s*b, 
            0, -δϕ*b, -δϕ*a, 
            -δθ*c, δϕ*c*a - δθ*s*b, -δϕ*b*c - δθ*a*s
        )

    elseif rot_seq == :YZ 
        return DCM(
            -a*s*δθ-b*c*δϕ, b*s*δθ-a*c*δϕ, c*δθ, 
            a*δϕ, -b*δϕ, 0, 
            b*s*δϕ-a*c*δθ, b*c*δθ+a*s*δϕ, -s*δθ
        )

    elseif rot_seq == :ZX 
        return DCM(
            -s*δθ, b*s*δϕ-a*c*δθ, b*c*δθ+a*s*δϕ, 
            c*δθ, -b*c*δϕ-a*s*δθ, b*s*δθ-a*c*δϕ, 
            0, a*δϕ, -b*δϕ        
        )

    elseif rot_seq == :ZY 
        return DCM(
            -a*s*δθ-b*c*δϕ, -c*δθ, a*c*δϕ-b*s*δθ, 
            a*c*δθ -b*s*δϕ, -s*δθ, b*c*δθ+a*s*δϕ, 
            -a*δϕ, 0, -b*δϕ
        )
    
    else 
        throw(ArgumentError(
            "The rotation sequence :$rot_seq is not valid."
        ))
    end

end

function angle_to_δdcm(θ, ϕ, γ, rot_seq::Symbol)

    s, c = sincos(θ[1])
    b, a = sincos(ϕ[1])
    e, d = sincos(γ[1])

    δγ, δϕ, δθ = γ[2], ϕ[2], θ[2]

    if rot_seq == :ZYX
        return DCM(
            -δθ*a*s         - δϕ*b*c, 
            -δθ*(d*c+e*b*s) + δϕ*e*a*c  + δγ*(e*s+d*b*c),
             δθ*(e*c-d*b*s) + δϕ*d*a*c  + δγ*(d*s-e*b*c),
             δθ*a*c         - δϕ*b*s,
             δθ*(e*b*c-d*s) + δϕ*e*a*s  + δγ*(d*b*s-e*c),
             δθ*(e*s+d*b*c) + δϕ*a*d*s  - δγ*(d*c+e*b*s),
                            - δϕ*a, 
                            - δϕ*e*b    + δγ*d*a, 
                            - δϕ*d*b    - δγ*e*a
         )

    elseif rot_seq == :XYX
        return DCM(
                              -δϕ*b, 
                              δϕ*e*a    + δγ*d*b, 
                              δϕ*d*a    - δγ*e*b, 
             δθ*b*c         + δϕ*a*s, 
            -δθ*(d*s+e*a*c) + δϕ*e*b*s  - δγ*(e*c+d*a*s),
             δθ*(e*s-d*a*c) + δϕ*d*b*s  + δγ*(e*a*s-d*c),
             δθ*b*s         - δϕ*a*c, 
             δθ*(d*c-e*a*s) - δϕ*e*b*c  + δγ*(d*a*c-e*s),
            -δθ*(e*c+d*a*s) - δϕ*d*b*c  - δγ*(d*s+e*a*c)
        )

    elseif rot_seq == :XYZ
        return DCM(
                            - δϕ*d*b    - δγ*e*a, 
                              δϕ*e*b    - δγ*d*a, 
                              δϕ*a, 
             δθ*(d*b*c-e*s) + δϕ*d*a*s  + δγ*(d*c-e*b*s),
            -δθ*(e*b*c+d*s) - δϕ*e*a*s  - δγ*(d*b*s+e*c),
            -δθ*a*c         + δϕ*b*s,
             δθ*(d*b*s+e*c) - δϕ*d*a*c  + δγ*(e*b*c+d*s),
             δθ*(d*c-e*b*s) + δϕ*e*a*c  + δγ*(d*b*c-e*s),
            -δθ*a*s         - δϕ*b*c
        )

    elseif rot_seq == :XZX
        return DCM(
                            - δϕ*b, 
                            - δϕ*d*a    + δγ*e*b, 
                              δϕ*e*a    + δγ*d*b, 
            -δθ*b*s         + δϕ*a*c,
            -δθ*(d*a*s+e*c) - δϕ*d*b*c  - δγ*(e*a*c+d*s),
             δθ*(e*a*s-d*c) + δϕ*e*b*c  + δγ*(e*s-d*a*c),
             δθ*b*c         + δϕ*a*s, 
             δθ*(d*a*c-e*s) - δϕ*d*b*s  + δγ*(d*c-e*a*s),
            -δθ*(e*a*c+d*s) + δϕ*e*b*s  - δγ*(d*a*s+e*c)
        )

    elseif rot_seq == :XZY
        return DCM(
                            - δϕ*d*b      - δγ*e*a, 
                            - δϕ*a,  
                            - δϕ*e*b    + δγ*d*a, 
             δθ*(e*c-d*b*s) + δϕ*d*a*c    + δγ*(d*s-e*b*c), 
            -δθ*a*s         - δϕ*b*c,  
            -δθ*(d*c+e*b*s) + δϕ*e*a*c  + δγ*(d*b*c+e*s),
             δθ*(d*b*c+e*s) + δϕ*d*a*s    - δγ*(e*b*s+d*c),
             δθ*a*c         - δϕ*b*s, 
             δθ*(e*b*c-d*s) + δϕ*e*a*s  + δγ*(d*b*s-e*c)
        )

    elseif rot_seq == :YXY
        return DCM(
            -δθ*(d*s+e*a*c) + δϕ*e*b*s  - δγ*(e*c+d*a*s),
             δθ*b*c         + δϕ*a*s,
             δθ*(d*a*c-e*s) - δϕ*d*b*s  + δγ*(d*c-e*a*s),
                              δϕ*e*a    + δγ*d*b,
                            - δϕ*b,
                            - δϕ*d*a    + δγ*e*b,
             δθ*(e*a*s-d*c) + δϕ*e*b*c  + δγ*(e*s-d*a*c),
            -δθ*b*s         + δϕ*a*c,
            -δθ*(e*c+d*a*s) - δϕ*d*b*c  - δγ*(d*s+e*a*c)

        )

    elseif rot_seq == :YXZ
        return DCM(
             δθ*(e*b*c-d*s) + δϕ*e*a*s  + δγ*(d*b*s-e*c),
             δθ*(e*s+d*b*c) + δϕ*d*a*s  - δγ*(d*c+e*b*s),
             δθ*a*c         - δϕ*b*s, 
                            - δϕ*e*b    + δγ*d*a,
                            - δϕ*d*b    - δγ*e*a,
                            - δϕ*a, 
            -δθ*(d*c+e*b*s) + δϕ*e*a*c  + δγ*(e*s+d*b*c),
             δθ*(e*c-d*b*s) + δϕ*d*a*c  + δγ*(d*s-e*b*c),
            -δθ*a*s         - δϕ*b*c, 
        )

    elseif rot_seq == :YZX
        return DCM(
            -δθ*a*s         - δϕ*b*c, 
             δθ*(d*b*s+e*c) - δϕ*d*a*c  + δγ*(e*b*c+d*s),
             δθ*(d*c-e*b*s) + δϕ*e*a*c  + δγ*(d*b*c-e*s),
                              δϕ*a,
                            - δϕ*d*b    - δγ*e*a,
                              δϕ*e*b    - δγ*d*a,
            -δθ*a*c         + δϕ*b*s,
             δθ*(d*b*c-e*s) + δϕ*d*a*s  + δγ*(d*c-e*b*s),
            -δθ*(e*b*c+d*s) - δϕ*e*a*s  - δγ*(d*b*s+e*c)
        )

    elseif rot_seq == :YZY
        return DCM(
            -δθ*(e*c+d*a*s) - δϕ*d*b*c  - δγ*(e*a*c+d*s),
             δθ*b*s         - δϕ*a*c, 
             δθ*(d*c-e*a*s) - δϕ*e*b*c  + δγ*(d*a*c-e*s),
                              δϕ*d*a    - δγ*e*b,
                            - δϕ*b, 
                              δϕ*e*a    + δγ*d*b,
             δθ*(e*s-d*a*c) + δϕ*d*b*s  + δγ*(e*a*s-d*c),
             δθ*b*c         + δϕ*a*s, 
            -δθ*(e*a*c+d*s) + δϕ*e*b*s  - δγ*(d*a*s+e*c)
        )

    elseif rot_seq == :ZXY
        return DCM(
            -δθ*(d*s+e*b*c) - δϕ*e*a*s  - δγ*(e*c+d*b*s),
            -δθ*a*c         + δϕ*b*s,
             δθ*(d*b*c-e*s) + δϕ*d*a*s  + δγ*(d*c-e*b*s),
             δθ*(d*c-e*b*s) + δϕ*e*a*c  + δγ*(d*b*c-e*s),
            -δθ*a*s         - δϕ*b*c,
             δθ*(e*c+d*b*s) - δϕ*d*a*c  + δγ*(d*s+e*b*c),
                              δϕ*e*b    - δγ*d*a,
                              δϕ*a,
                            - δϕ*d*b    - δγ*e*a
        )

    elseif rot_seq == :ZXZ
        return DCM(
            -δθ*(d*s+e*a*c) + δϕ*e*b*s  - δγ*(e*c+d*a*s),
             δθ*(e*s-d*a*c) + δϕ*d*b*s  + δγ*(e*a*s-d*c),
             δθ*b*c         + δϕ*a*s,
             δθ*(d*c-e*a*s) - δϕ*e*b*c  + δγ*(d*a*c-e*s),
            -δθ*(e*c+d*a*s) - δϕ*d*b*c  - δγ*(d*s+e*a*c),
             δθ*b*s         - δϕ*a*c,
                              δϕ*a*e    + δγ*d*b,
                              δϕ*d*a    - δγ*e*b,
                            - δϕ*b
        )
        
    elseif rot_seq == :ZYZ
        return DCM(
            -δθ*(d*a*s+e*c) - δϕ*d*b*c  - δγ*(e*a*c+d*s),
             δθ*(e*a*s-d*c) + δϕ*e*b*c  + δγ*(e*s-d*a*c),
            -δθ*b*s         + δϕ*a*c,
             δθ*(d*a*c-e*s) - δϕ*d*b*s  + δγ*(d*c-e*a*s),
            -δθ*(e*a*c+d*s) + δϕ*e*b*s  - δγ*(d*a*s+e*c),
             δθ*b*c         + δϕ*a*s,
                            - δϕ*a*d    + δγ*e*b,
                              δϕ*e*a    + δγ*d*b,
                            - δϕ*b
        )

    else
        throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
    end
end
