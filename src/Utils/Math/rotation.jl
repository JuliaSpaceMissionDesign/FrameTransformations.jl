
export skew 

"""
    skew(a)

Create a skew matrix from the first three elements of `a`.
"""
skew(a) = SMatrix{3, 3}(0, a[3], -a[2], -a[3], 0, a[1], a[2], -a[1], 0)


"""
    angle_to_δdcm(θ[, ϕ[, γ]], rot_seq::Symbol = :ZYX)

Compute the derivative of the direction cosine matrix that perform a set of rotations (`θ`, 
`ϕ`, `γ`) about the coordinate axes in `rot_seq`. Each rotation input must be an indexable 
objected which includes the angle and its first time derivative.

The rotation sequence is defined by a `Symbol` specifing the rotation axes. The possible 
values depend on the number of rotations as follows: 

- **1 rotation** (`θ₁`): `:X`, `:Y`, or `:Z`.
- **2 rotations** (`θ₁`, `θ₂`): `:XY`, `:XZ`, `:YX`, `:YZ`, `:ZX`, or `:ZY`.
- **3 rotations** (`θ₁`, `θ₂`, `θ₃`): `:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`,
    `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`, `:ZYX`, or `:ZYZ`

!!! note 
    This function assigns `dcm = A3 * A2 * A1` in which `Ai` is the DCM related with
    the *i*-th rotation, `i Є [1,2,3]`. If the *i*-th rotation is not specified,
    then `Ai = I`.

### See also 
See also [`angle_to_δ²dcm`](@ref)
"""
function angle_to_δdcm(θ, rot_seq::Symbol=:ZYX)
    
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
    @inbounds _3angles_to_δdcm((θ[1], ϕ[1], γ[1], θ[2], ϕ[2], γ[2]), rot_seq)
end

"""
    _3angles_to_δdcm(θ, rot_seq::Symbol)

Compute the time derivative of the DCM with all the angles stored in a single vector 
to optimise computations. 
"""
function _3angles_to_δdcm(θ, rot_seq::Symbol)

    s, c = sincos(θ[1])
    b, a = sincos(θ[2])
    e, d = sincos(θ[3])

    δθ, δϕ, δγ = θ[4], θ[5], θ[6]

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
                            - δϕ*d*b    - δγ*e*a, 
                            - δϕ*a,  
                            - δϕ*e*b    + δγ*d*a, 
             δθ*(e*c-d*b*s) + δϕ*d*a*c  + δγ*(d*s-e*b*c), 
            -δθ*a*s         - δϕ*b*c,  
            -δθ*(d*c+e*b*s) + δϕ*e*a*c  + δγ*(d*b*c+e*s),
             δθ*(d*b*c+e*s) + δϕ*d*a*s  - δγ*(e*b*s+d*c),
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


""" 
    angle_to_δ²dcm(θ[, ϕ[, γ]], rot_seq::Symbol = :ZYX)

Compute the second order time derivative of the direction cosine matrix that perform 
a set of rotations (`θ`, `ϕ`, `γ`) about the coordinate axes in `rot_seq`. Each rotation 
input must be an indexable objected which includes the angle and its first and second order 
time derivatives.

The rotation sequence is defined by a `Symbol` specifing the rotation axes. The possible 
values depend on the number of rotations as follows: 

- **1 rotation** (`θ₁`): `:X`, `:Y`, or `:Z`.
- **2 rotations** (`θ₁`, `θ₂`): `:XY`, `:XZ`, `:YX`, `:YZ`, `:ZX`, or `:ZY`.
- **3 rotations** (`θ₁`, `θ₂`, `θ₃`): `:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`,
    `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`, `:ZYX`, or `:ZYZ`

!!! note 
    This function assigns `dcm = A3 * A2 * A1` in which `Ai` is the DCM related with
    the *i*-th rotation, `i Є [1,2,3]`. If the *i*-th rotation is not specified,
    then `Ai = I`.

### See also 
See also [`angle_to_δdcm`](@ref)
"""
function angle_to_δ²dcm(θ, rot_seq::Symbol=:ZYX)

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


function angle_to_δ²dcm(θ, ϕ, rot_seq::Symbol)

    s, c = sincos(θ[1])
    b, a = sincos(ϕ[1])

    δϕ, δ²ϕ = ϕ[2], ϕ[3]
    δθ, δ²θ = θ[2], θ[3]

    δθ², δϕ² = δθ^2, δϕ^2

    δϕθ = 2*δϕ*δθ
    B = (δθ²+ δϕ²)

    if rot_seq == :XY
        return DCM(
                        - δ²ϕ*b     - δϕ²*a, 
            0,
                          δ²ϕ*a     - δϕ²*b,
             δ²θ*c*b    + δ²ϕ*a*s   - B*b*s     + δϕθ*a*c,
            -δ²θ*s                  - δθ²*c,
            -δ²θ*a*c    + δ²ϕ*b*s   + B*a*s     + δϕθ*c*b,
             δ²θ*b*s    - δ²ϕ*a*c   + B*c*b     + δϕθ*a*s,
             δ²θ*c                  - δθ²*s,
            -δ²θ*a*s    - δ²ϕ*c*b   - B*a*c     + δϕθ*b*s
        )

    elseif rot_seq == :XZ
        return DCM(
                        - δ²ϕ*b     - δϕ²*a,
                        - δ²ϕ*a     + δϕ²*b,
            0,
            -δ²θ*b*s    + δ²ϕ*a*c   - B*c*b     - δϕθ*a*s,
            -δ²θ*a*s    - δ²ϕ*c*b   - B*a*c     + δϕθ*b*s,
            -δ²θ*c                  + δθ²*s ,
             δ²θ*c*b    + δ²ϕ*a*s   - B*b*s     + δϕθ*a*c,
             δ²θ*a*c    - δ²ϕ*b*s   - B*a*s     - δϕθ*c*b,
            -δ²θ*s                  - δθ²*c
        )
        
        
    elseif rot_seq == :YX 
        return DCM(
            -δ²θ*s                  - δθ²*c,
             δ²θ*c*b    + δ²ϕ*a*s   - B*b*s     + δϕθ*a*c,
             δ²θ*a*c    - δ²ϕ*b*s   - B*a*s     - δϕθ*c*b,
            0,
                        - δ²ϕ*b     - δϕ²*a,
                        - δ²ϕ*a     + δϕ²*b,
            -δ²θ*c                  + δθ²*s,
            -δ²θ*b*s    + δ²ϕ*a*c   - B*c*b     - δϕθ*a*s,
            -δ²θ*a*s    - δ²ϕ*c*b   - B*a*c     + δϕθ*b*s 
        )

    elseif rot_seq == :YZ 
        return DCM(
            -δ²θ*a*s    - δ²ϕ*c*b   - B*a*c     + δϕθ*b*s,
             δ²θ*b*s    - δ²ϕ*a*c   + B*c*b     + δϕθ*a*s,
             δ²θ*c                  - δθ²*s,
                          δ²ϕ*a     - δϕ²*b,
                        - δ²ϕ*b     - δϕ²*a ,
            0,
            -δ²θ*a*c    + δ²ϕ*b*s   + B*a*s     + δϕθ*c*b,
            +δ²θ*c*b    + δ²ϕ*a*s   - B*b*s     + δϕθ*a*c,
            -δ²θ*s                  - δθ²*c
        )

    elseif rot_seq == :ZX 
        return DCM(
            -δ²θ*s                  - δθ²*c,
            -δ²θ*a*c    + δ²ϕ*b*s   + B*a*s     + δϕθ*c*b,
             δ²θ*c*b    + δ²ϕ*a*s   - B*b*s     + δϕθ*a*c,
             δ²θ*c                  - δθ²*s,
            -δ²θ*a*s    - δ²ϕ*c*b   - B*a*c     + δϕθ*b*s,
             δ²θ*b*s    - δ²ϕ*a*c   + B*c*b     + δϕθ*a*s,
            0,
                          δ²ϕ*a     - δϕ²*b,
                        - δ²ϕ*b     - δϕ²*a
        )

    elseif rot_seq == :ZY 
        return DCM(
            -δ²θ*a*s    - δ²ϕ*c*b   - B*a*c     + δϕθ*b*s,
            -δ²θ*c                  + δθ²*s,
            -δ²θ*b*s    + δ²ϕ*a*c   - B*c*b     - δϕθ*a*s,
             δ²θ*a*c    - δ²ϕ*b*s   - B*a*s     - δϕθ*c*b,
            -δ²θ*s                  - δθ²*c,
             δ²θ*c*b    + δ²ϕ*a*s   - B*b*s + δϕθ*a*c,
                        - δ²ϕ*a     + δϕ²*b,
            0,
                        - δ²ϕ*b     - δϕ²*a 
        )
    
    else 
        throw(ArgumentError(
            "The rotation sequence :$rot_seq is not valid."
        ))
    end

end


function angle_to_δ²dcm(θ, ϕ, γ, rot_seq::Symbol)
    @inbounds _3angles_to_δ²dcm((θ[1], ϕ[1], γ[1], 
                     θ[2], ϕ[2], γ[2], 
                     θ[3], ϕ[3], γ[3]), rot_seq)
end


"""
    _3angles_to_δ²dcm(θ, rot_seq::Symbol)

Compute the 2nd time derivative of the DCM with all the angles stored in a single vector 
to optimise computations. 
"""
function _3angles_to_δ²dcm(θ, rot_seq::Symbol)

    s, c = sincos(θ[1])
    b, a = sincos(θ[2])
    e, d = sincos(θ[3])

    δθ, δϕ, δγ    = θ[4], θ[5], θ[6]
    δ²θ, δ²ϕ, δ²γ = θ[7], θ[8], θ[9]

    δθ², δϕ², δγ² = δθ^2, δϕ^2, δγ^2

    δϕθ = 2*δϕ*δθ
    δϕγ = 2*δϕ*δγ
    δθγ = 2*δθ*δγ

    A = (δθ² + δϕ² + δγ²)
    B = δθ² + δϕ²
    C = δθ² + δγ²
    D = δγ² + δϕ²

	if rot_seq == :ZYX
		return DCM(
			-a*(c*B + s*δ²θ) - b*(c*δ²ϕ - s*δϕθ),
            c*(d*(a*δϕγ + δ²γ*b - δ²θ) + e*(a*δ²ϕ - A*b + δθγ)) + s*(e*(δ²γ - a*δϕθ - δ²θ*b) + d*(C-δθγ*b)),
            c*(a*(d*δ²ϕ - e*δϕγ) + e*(δ²θ - δ²γ*b) + d*(δθγ- A*b)) + s*(d*(δ²γ-a*δϕθ -δ²θ*b) - e*(C-δθγ*b)),
			a*(c*δ²θ - B*s) - b*(c*δϕθ + s*δ²ϕ),
            c*(e*(a*δϕθ -δ²γ + b*δ²θ) + d*(δθγ*b-C)) + s*(d*(a*δϕγ + δ²γ*b - δ²θ) +e*(a*δ²ϕ + δθγ - A*b)),
            c*d*(a*δϕθ-δ²γ+δ²θ*b) + d*s*(a*δ²ϕ- A*b + δθγ) -e*s*(a*δϕγ - δ²θ + δ²γ*b) + c*e*(C-δθγ*b),
			-a*δ²ϕ + b*δϕ²,
            a*(d*δ²γ-D*e) - b*(d*δϕγ + e*δ²ϕ),
			-a*(D*d + e*δ²γ) + b*(-d*δ²ϕ + e*δϕγ),
		)
	elseif rot_seq == :XYX
		return DCM(
			-a*δϕ² - b*δ²ϕ,
			a*(d*δϕγ + e*δ²ϕ) + b*(d*δ²γ - D*e),
			a*(d*δ²ϕ - e*δϕγ) - b*(D*d + e*δ²γ),
			a*(c*δϕθ + s*δ²ϕ) + b*(c*δ²θ - B*s),
            e*(s*(A*a + b*δ²ϕ + δθγ) + c*(b*δϕθ - δ²γ - δ²θ*a)) + d*(s*(b*δϕγ - δ²γ*a - δ²θ) -c*(C + δθγ*a)),
            d*(s*(A*a+b*δ²ϕ+δθγ) + c*(b*δϕθ-δ²γ-a*δ²θ)) + e*(s*(δ²γ*a- b*δϕγ+δ²θ) + c*(C+δθγ*a)),
			a*(s*δϕθ-c*δ²ϕ ) + b*(B*c + s*δ²θ),
            c*(e*(-A*a-b*δ²ϕ-δθγ) + d*(δ²θ+δ²γ*a- b*δϕγ)) + s*(e*(b*δϕθ-δ²γ-a*δ²θ) - d*(C+a*δθγ)),
            c*(d*(-A*a-δθγ) + b*(-d*δ²ϕ+e*δϕγ ) - e*(δ²θ+δ²γ*a)) + s*(d*(b*δϕθ-δ²γ- δ²θ*a) + e*(C+δθγ*a)),
		)
	elseif rot_seq == :XYZ 
		return DCM(
            e*(b*δϕγ-a*δ²γ) - d*(D*a+b*δ²ϕ),
            a*(e*D-d*δ²γ) + b*(d*δϕγ + e*δ²ϕ),
			a*δ²ϕ - b*δϕ²,
            d*(c*(a*δϕθ+δ²γ+b*δ²θ) + s*(a*δ²ϕ-A*b-δθγ)) + e*(s*(-a*δϕγ-b*δ²γ-δ²θ) - c*(C+δθγ*b)),
            - c*(d*(C+b*δθγ)+e*(δϕθ*a+δ²γ+b*δ²θ)) + s*(-d*(a*δϕγ+b*δ²γ+δ²θ) + e*(δθγ+A*b-a*δ²ϕ)),
			-a*(c*δ²θ - s*B) + b*(c*δϕθ + s*δ²ϕ),
            c*(a*(-d*δ²ϕ+e*δϕγ) + e*(δ²γ*b+δ²θ) + d*(A*b+δθγ)) + s*(-e*(C+b*δθγ) + d*(a*δϕθ+δ²γ+b*δ²θ)),
            c*(a*(d*δϕγ) + e*(δ²ϕ*a-A*b-δθγ) + d*(δ²γ*b+δ²θ)) - s*(e*(δϕθ*a+δ²γ+b*δ²θ) + d*(C+b*δθγ)),
			-a*(c*B + s*δ²θ) - b*(c*δ²ϕ - s*δϕθ),
		)
	elseif rot_seq == :XZX 
		return DCM(
			-a*δϕ² - b*δ²ϕ,
			a*(-d*δ²ϕ + e*δϕγ) + b*(D*d + e*δ²γ),
			a*(d*δϕγ + e*δ²ϕ) + b*(d*δ²γ - D*e),
			a*(c*δ²ϕ - s*δϕθ) - b*(B*c + s*δ²θ),
            c*(e*(b*δϕγ-a*δ²γ-δ²θ) - d*(A*a+b*δ²ϕ+δθγ)) + s*(d*(b*δϕθ-δ²γ-a*δ²θ) + e*(C+a*δθγ)),
            c*(d*(b*δϕγ-a*δ²γ-δ²θ) + e*(A*a+b*δ²ϕ+δθγ)) + s*(e*(δ²γ-b*δϕθ+a*δ²θ) + d*(C+a*δθγ)),
			a*(c*δϕθ + s*δ²ϕ) + b*(c*δ²θ - B*s),
            d*(-s*(+A*a+b*δ²ϕ+δθγ) + c*(a*δ²θ+δ²γ-b*δϕθ)) + e*(s*(b*δϕγ-a*δ²γ-δ²θ) - c*(C+a*δθγ)),
            e*(s*(A*a+b*δ²ϕ+δθγ) + c*(b*δϕθ-δ²γ-a*δ²θ)) + d*(s*(b*δϕγ-a*δ²γ-δ²θ) - c*(C+a*δθγ)),
		)
	elseif rot_seq == :XZY 
		return DCM(
			-a*(d*D + e*δ²γ) + b*(-d*δ²ϕ + e*δϕγ),
			-a*δ²ϕ + b*δϕ²,
			a*(d*δ²γ - e*D) - b*(d*δϕγ + e*δ²ϕ),
            c*(d*(δ²ϕ*a-b*A+δθγ) + e*(δ²θ-a*δϕγ-b*δ²γ)) + s*(e*(δθγ*b-C) + d*(δ²γ-a*δϕθ-b*δ²θ)),
			-a*(c*B + s*δ²θ) + b*(-c*δ²ϕ + s*δϕθ),
            c*(d*(a*δϕγ+b*δ²γ-δ²θ) + e*(a*δ²ϕ-A*b+δθγ)) + s*(e*(δ²γ-a*δϕθ-b*δ²θ) + d*(C-b*δθγ)),
            c*(d*(a*δϕθ-δ²γ+b*δ²θ) + e*(C-b*δθγ)) + s*(e*(δ²θ-a*δϕγ-b*δ²γ) + d*(δθγ+a*δ²ϕ-b*A)),
			a*(c*δ²θ - s*B) - b*(c*δϕθ + s*δ²ϕ),
            c*(d*(δθγ*b-C) + e*(a*δϕθ-δ²γ+b*δ²θ)) + s*(d*(a*δϕγ+b*δ²γ-δ²θ)+e*(a*δ²ϕ-b*A+δθγ)),
		)

	elseif rot_seq == :YXY 
		return DCM(
            c*(-d*(C+a*δθγ)+e*(b*δϕθ-δ²γ-a*δ²θ)) + s*(d*(b*δϕγ-a*δ²γ-δ²θ)+e*(A*a+b*δ²ϕ+δθγ)),
			a*(c*δϕθ + s*δ²ϕ) + b*(c*δ²θ - s*B),
            c*(d*(a*δ²θ+δ²γ-b*δϕθ) - e*(C+a*δθγ)) + s*(-d*(a*A+b*δ²ϕ+δθγ) + e*(-δ²θ-a*δ²γ+b*δϕγ)),
			a*(d*δϕγ + e*δ²ϕ) + b*(d*δ²γ - e*D),
			-a*δϕ² - b*δ²ϕ,
			a*(-d*δ²ϕ + e*δϕγ) + b*(d*D + e*δ²γ),
            c*(d*(b*δϕγ-a*δ²γ-δ²θ)+e*(δθγ+a*A+b*δ²ϕ)) + s*(d*(C+a*δθγ)+e*(a*δ²θ+δ²γ-b*δϕθ)),
            a*(c*δ²ϕ - s*δϕθ) - b*(c*B + s*δ²θ),
            c*(-d*(+a*A+b*δ²ϕ+δθγ)+e*(b*δϕγ-a*δ²γ-δ²θ)) + s*(d*(b*δϕθ-δ²γ-a*δ²θ) + e*(C+a*δθγ)),
		)

	elseif rot_seq == :YXZ 
		return DCM(
            c*(d*(b*δθγ-C) + e*(a*δϕθ-δ²γ+b*δ²θ)) + s*(d*(a*δϕγ+b*δ²γ-δ²θ) + e*(a*δ²ϕ-b*A+δθγ)),
            c*(d*(a*δϕθ-δ²γ+b*δ²θ) + e*(C-b*δθγ)) + s*(d*(a*δ²ϕ-b*A+δθγ) + e*(δ²θ-a*δϕγ-b*δ²γ)),
			a*(c*δ²θ - s*B) - b*(c*δϕθ + s*δ²ϕ),
			a*(d*δ²γ - e*D) - b*(d*δϕγ + e*δ²ϕ),
			-a*(d*D + e*δ²γ) + b*(-d*δ²ϕ + e*δϕγ),
			-a*δ²ϕ + b*δϕ²,
            c*(d*(a*δϕγ+b*δ²γ-δ²θ) + e*(a*δ²ϕ-b*A+δθγ)) + s*(d*(C-b*δθγ) + e*(δ²γ-a*δϕθ-b*δ²θ)),
            c*(d*(a*δ²ϕ-b*A+δθγ) + e*(δ²θ-a*δϕγ-b*δ²γ)) + s*(d*(δ²γ-a*δϕθ-b*δ²θ) + e*(δθγ*b-C)),
			-a*(c*B + s*δ²θ) - b*(c*δ²ϕ - s*δϕθ),
		)
	elseif rot_seq == :YZX 
		return DCM(
			-a*(c*B +s*δ²θ) + b*(-c*δ²ϕ + s*δϕθ),
            c*(d*(-a*δ²ϕ+b*A+δθγ) + e*(δ²θ+a*δϕγ+b*δ²γ)) + s*(d*(a*δϕθ+δ²γ+b*δ²θ) + e*(-C-b*δθγ)),
            c*(d*(a*δϕγ+b*δ²γ+δ²θ) + e*(a*δ²ϕ-b*A-δθγ)) + s*(-d*(C+b*δθγ) + e*(-b*δ²θ-δ²γ-a*δϕθ)),
			a*δ²ϕ - b*δϕ²,
			-a*(d*D + e*δ²γ) + b*(-d*δ²ϕ + e*δϕγ),
			a*(-d*δ²γ + e*D) + b*(d*δϕγ + e*δ²ϕ),
			a*(-c*δ²θ + s*B) + b*(c*δϕθ + s*δ²ϕ),
            c*(d*(a*δϕθ+δ²γ+b*δ²θ) - e*(C+b*δθγ)) + s*(d*(a*δ²ϕ-b*A-δθγ) - e*(δ²θ+a*δϕγ+b*δ²γ)),
            c*(-d*(C+b*δθγ) - e*(a*δϕθ+δ²γ+b*δ²θ)) + s*(-d*(δ²θ+a*δϕγ+b*δ²γ) + e*(δθγ-a*δ²ϕ+b*A)),
		)
	elseif rot_seq == :YZY 
		return DCM(
            c*(d*(-a*A-b*δ²ϕ-δθγ) + e*(b*δϕγ-a*δ²γ-δ²θ)) + s*(d*(b*δϕθ-δ²γ-a*δ²θ) + e*(C+a*δθγ)),
			a*(-c*δ²ϕ + s*δϕθ) + b*(c*B + s*δ²θ),
            c*(d*(a*δ²γ-b*δϕγ+δ²θ) -e*(δθγ+a*A+b*δ²ϕ)) + s*(-d*(C+a*δθγ) + e*(b*δϕθ-δ²γ-a*δ²θ)),
			a*(d*δ²ϕ - e*δϕγ) - b*(d*D + e*δ²γ),
			-a*δϕ² - b*δ²ϕ,
			a*(d*δϕγ + e*δ²ϕ) + b*(d*δ²γ - e*D),
            c*(d*(b*δϕθ-δ²γ-a*δ²θ) + e*(C+a*δθγ)) + s*(d*(a*A+b*δ²ϕ+δθγ) + e*(δ²θ+δ²γ*a-b*δϕγ)),
			a*(c*δϕθ + s*δ²ϕ) + b*(c*δ²θ - s*B),
            c*(-d*(C+a*δθγ) + e*(b*δϕθ-δ²γ-a*δ²θ)) + s*(d*(b*δϕγ-a*δ²γ-δ²θ) + e*(a*A+b*δ²ϕ+δθγ)),
		)
	elseif rot_seq == :ZXY 
		return DCM(
            c*(-d*(C+b*δθγ) - e*(a*δϕθ+δ²γ+b*δ²θ)) + s*(-d*(δ²θ+b*δ²γ+a*δϕγ) + e*(δθγ+b*A-a*δ²ϕ)),
			a*(-c*δ²θ +s*B) + b*(c*δϕθ + s*δ²ϕ),
            c*(d*(a*δϕθ+δ²γ+b*δ²θ) - e*(C+b*δθγ)) + s*(d*(a*δ²ϕ-b*A-δθγ) - e*(b*δ²γ+a*δϕγ+δ²θ)),
            c*(d*(a*δϕγ+b*δ²γ+δ²θ) + e*(a*δ²ϕ-b*A-δθγ)) + s*(-d*(C+b*δθγ) - e*(a*δϕθ+δ²γ+b*δ²θ)),
			-a*(c*B +s*δ²θ) + b*(-c*δ²ϕ + s*δϕθ),
            c*(d*(-a*δ²ϕ+b*A+δθγ) + e*(a*δϕγ+b*δ²γ+δ²θ)) + s*(d*(a*δϕθ+δ²γ+b*δ²θ) - e*(C+b*δθγ)),
			a*(-d*δ²γ + e*D) + b*(d*δϕγ + e*δ²ϕ),
			a*δ²ϕ - b*δϕ²,
			-a*(d*D+e*δ²γ) + b*(-d*δ²ϕ + e*δϕγ),
		)
	elseif rot_seq == :ZXZ 
		return DCM(
            c*(-d*(C+a*δθγ) + e*(b*δϕθ-δ²γ-a*δ²θ)) + s*(d*(b*δϕγ-a*δ²γ-δ²θ) + e*(δθγ+a*A+b*δ²ϕ)),
            c*(d*(b*δϕθ-δ²γ-a*δ²θ) + e*(C+a*δθγ)) + s*(d*(a*A+b*δ²ϕ+δθγ) + e*(δ²θ-b*δϕγ+a*δ²γ)),
			a*(c*δϕθ + s*δ²ϕ) + b*(c*δ²θ - s*B),
            c*(d*(a*δ²γ-b*δϕγ+δ²θ) - e*(b*δ²ϕ+a*A+δθγ)) + s*(-d*(C+a*δθγ) + e*(b*δϕθ-δ²γ-a*δ²θ)),
            c*(-d*(a*A+b*δ²ϕ+δθγ) + e*(b*δϕγ-a*δ²γ-δ²θ)) + s*(d*(b*δϕθ-δ²γ-a*δ²θ) + e*(C+a*δθγ)),
			a*(-c*δ²ϕ + s*δϕθ) + b*(c*B + s*δ²θ),
			a*(d*δϕγ + e*δ²ϕ) + b*(d*δ²γ - e*D),
			a*(d*δ²ϕ - e*δϕγ) - b*(d*D + e*δ²γ),
			-a*δϕ² - b*δ²ϕ,
		)
	elseif rot_seq == :ZYZ 
		return DCM(
            c*(-d*(δθγ+a*A+b*δ²ϕ) + e*(b*δϕγ-a*δ²γ-δ²θ)) + s*(d*(b*δϕθ-δ²γ-a*δ²θ) + e*(C+a*δθγ)),
            c*(d*(b*δϕγ-a*δ²γ-δ²θ) + e*(a*A+b*δ²ϕ+δθγ)) + s*(d*(C+a*δθγ) + e*(a*δ²θ+δ²γ-b*δϕθ)),
			a*(c*δ²ϕ - s*δϕθ) - b*(c*B +s*δ²θ),
            c*(d*(a*δ²θ+δ²γ-b*δϕθ) - e*(C+a*δθγ)) + s*(-d*(δθγ+a*A+b*δ²ϕ) + e*(b*δϕγ-a*δ²γ-δ²θ)),
            c*(-d*(C+a*δθγ) + e*(b*δϕθ-δ²γ-a*δ²θ)) + s*(d*(b*δϕγ-a*δ²γ-δ²θ) + e*(a*A+b*δ²ϕ+δθγ)),
			a*(c*δϕθ + s*δ²ϕ) + b*(c*δ²θ - s*B),
			a*(-d*δ²ϕ + e*δϕγ) + b*(d*D + e*δ²γ),
			a*(d*δϕγ + e*δ²ϕ) + b*(d*δ²γ - e*D),
			-a*δϕ² - b*δ²ϕ,
		)
	else
		throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
	end


end


"""
    angle_to_δ³dcm(θ::Number[, ϕ::Number[, γ::Number]], rot_seq::Symbol = :Z)

Compute the second order time derivative of the direction cosine matrix that perform 
a set of rotations (`θ`, `ϕ`, `γ`) about the coordinate axes in `rot_seq`. Each rotation 
input must be an indexable objected which includes the angle and its time derivatives up to 
order 3 (jerk).

The rotation sequence is defined by a `Symbol` specifing the rotation axes. The possible 
values depend on the number of rotations as follows: 

- **1 rotation** (`θ₁`): `:X`, `:Y`, or `:Z`.
- **2 rotations** (`θ₁`, `θ₂`): `:XY`, `:XZ`, `:YX`, `:YZ`, `:ZX`, or `:ZY`.
- **3 rotations** (`θ₁`, `θ₂`, `θ₃`): `:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`,
    `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`, `:ZYX`, or `:ZYZ`

!!! note 
    This function assigns `dcm = A3 * A2 * A1` in which `Ai` is the DCM related with
    the *i*-th rotation, `i Є [1,2,3]`. If the *i*-th rotation is not specified,
    then `Ai = I`.

### See also 
See also [`angle_to_δdcm`](@ref)
"""
function angle_to_δ³dcm(θ, rot_seq::Symbol=:Z)

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

function angle_to_δ³dcm(θ, ϕ, rot_seq::Symbol)

	s, c = sincos(θ[1])
	b, a = sincos(ϕ[1])

	δθ, δ²θ, δ³θ = θ[2], θ[3], θ[4]
	δϕ, δ²ϕ, δ³ϕ = ϕ[2], ϕ[3], ϕ[4]

    δθ², δϕ² = δθ^2, δϕ^2
    δθ³, δϕ³ = δθ²*δθ, δϕ²*δϕ

	δθδ²ϕ = 3δθ*δ²ϕ 
	δθδ²θ = 3δθ*δ²θ 
	δϕδ²θ = 3δϕ*δ²θ 
	δϕδ²ϕ = 3δϕ*δ²ϕ 
	δϕδθ² = 3δϕ*δθ² 
	δϕ²δθ = 3δϕ²*δθ 

	A = δθδ²θ + δϕδ²ϕ
	B = δθδ²ϕ + δϕδ²θ
	C = -δ³ϕ + δϕ³ + δϕδθ²
	D = -δ³θ + δθ³ + δϕ²δθ
	E = δ³ϕ - δϕ³
	F = δ³θ - δθ³
		
	if rot_seq == :XY
		return DCM(
			-a*δϕδ²ϕ - b*E,0,
			a*E - b*δϕδ²ϕ,
			a*(c*B - s*C) - b*(c*D + s*A),
			-c*δθδ²θ - s*F,
			a*(c*D + s*A) + b*(c*B - s*C),
			a*(c*C + s*B) + b*(c*A - s*D),
			c*F - s*δθδ²θ,
			-a*(c*A - s*D) + b*(c*C + s*B),
		)
	elseif rot_seq == :XZ
		return DCM(
			-a*δϕδ²ϕ - b*E,
			-a*E + b*δϕδ²ϕ,
			0,
			-a*(c*C + s*B) - b*(c*A - s*D),
			-a*(c*A - s*D) + b*(c*C + s*B),
			-c*F + s*δθδ²θ,
			a*(c*B - s*C) - b*(c*D + s*A),
			-a*(c*D + s*A) - b*(c*B - s*C),
			-c*δθδ²θ - s*F,
		)
	elseif rot_seq == :YX
		return DCM(
			-c*δθδ²θ - s*F,
			a*(c*B - s*C) - b*(c*D + s*A),
			-a*(c*D + s*A) - b*(c*B - s*C),
			0,
			-a*δϕδ²ϕ - b*E,
			-a*E + b*δϕδ²ϕ,
			-c*F + s*δθδ²θ,
			-a*(c*C + s*B) - b*(c*A - s*D),
			-a*(c*A - s*D) + b*(c*C + s*B),
		)
	elseif rot_seq == :YZ
		return DCM(
			-a*(c*A - s*D) + b*(c*C + s*B),
			a*(c*C + s*B) + b*(c*A - s*D),
			c*F - s*δθδ²θ,
			a*E - b*δϕδ²ϕ,
			-a*δϕδ²ϕ - b*E,
			0,
			a*(c*D + s*A) + b*(c*B - s*C),
			a*(c*B - s*C) - b*(c*D + s*A),
			-c*δθδ²θ - s*F,
		)
	elseif rot_seq == :ZX
		return DCM(
			-c*δθδ²θ - s*F,
			a*(c*D + s*A) + b*(c*B - s*C),
			a*(c*B - s*C) - b*(c*D + s*A),
			c*F - s*δθδ²θ,
			-a*(c*A - s*D) + b*(c*C + s*B),
			a*(c*C + s*B) + b*(c*A - s*D),
			0,
			a*E - b*δϕδ²ϕ,
			-a*δϕδ²ϕ - b*E,
		)
	elseif rot_seq == :ZY
		return DCM(
			-a*(c*A - s*D) + b*(c*C + s*B),
			-c*F + s*δθδ²θ,
			-a*(c*C + s*B) - b*(c*A - s*D),
			-a*(c*D + s*A) - b*(c*B - s*C),
			-c*δθδ²θ - s*F,
			a*(c*B - s*C) - b*(c*D + s*A),
			-a*E + b*δϕδ²ϕ,
			0,
			-a*δϕδ²ϕ - b*E,
		)
	else
		throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
	end

end

function angle_to_δ³dcm(θ, ϕ, γ, rot_seq::Symbol)
    @inbounds _3angles_to_δ³dcm(
        (
            θ[1], ϕ[1], γ[1], 
            θ[2], ϕ[2], γ[2], 
            θ[3], ϕ[3], γ[3], 
            θ[4], ϕ[4], γ[4]
        ), 
        rot_seq
    )

end


function _3angles_to_δ³dcm(θ, rot_seq::Symbol)

    s, c = sincos(θ[1])
    b, a = sincos(θ[2])
    e, d = sincos(θ[3])

    δθ, δϕ, δγ    = θ[4], θ[5], θ[6]
    δ²θ, δ²ϕ, δ²γ = θ[7], θ[8], θ[9]
	δ³θ, δ³ϕ, δ³γ = θ[10], θ[11], θ[12]

    δθ², δϕ², δγ² = δθ^2, δϕ^2, δγ^2
    δθ³, δϕ³, δγ³ = δθ²*δθ, δϕ²*δϕ, δγ²*δγ

	δθδ²θ = 3δθ*δ²θ
	δθδ²ϕ = 3δθ*δ²ϕ
	δθδ²γ = 3δθ*δ²γ

	δγδ²θ = 3δγ*δ²θ
	δγδ²ϕ = 3δγ*δ²ϕ 
	δγδ²γ = 3δγ*δ²γ

	δϕδ²θ = 3δϕ*δ²θ
	δϕδ²ϕ = 3δϕ*δ²ϕ
	δϕδ²γ = 3δϕ*δ²γ

	δγδθ² = 3δγ*δθ²
	δγδϕ² = 3δγ*δϕ²
	δγ²δθ = 3δγ²*δθ
	δϕ²δθ = 3δϕ²*δθ
	δγ²δϕ = 3δγ²*δϕ
	δϕδθ² = 3δϕ*δθ²

	δγδϕδθ = 6*δγ*δϕ*δθ

	A = δθδ²θ + δϕδ²ϕ
	Z = δγδ²γ + A
	B = δγδ²γ + δϕδ²ϕ
	C = δγδ²γ + δθδ²θ
	D = δγδ²ϕ + δϕδ²γ
	E = δθδ²ϕ + δϕδ²θ
	F = δγδ²θ + δθδ²γ

	G = -δ³θ + δθ³ + δϕ²δθ
	N = -δ³θ + δθ³ + δγ²δθ 
	H = -δ³ϕ + δϕ³ + δϕδθ²
	I = -δ³ϕ + δϕ³ + δγ²δϕ 
	L = -δ³γ + δγ³ + δγδϕ²
	M = -δ³γ + δγ³ + δγδθ²

	if rot_seq == :ZYX
		return DCM(
			-a*(c*A - s*G) + b*(c*H + s*E),
			a*(c*(d*D - e*(I + δϕδθ²)) - s*(d*δγδϕδθ + e*E)) - b*(c*(d*(M + δγδϕ²) + e*Z) + s*(d*F - e*(N + δϕ²δθ))) + c*(d*N + e*F) + s*(d*C - e*M),
			-a*(c*(d*(I + δϕδθ²) + e*D) + s*(d*E - e*δγδϕδθ)) - b*(c*(d*Z - e*(M + δγδϕ²)) - s*(d*(N + δϕ²δθ) + e*F)) + c*(d*F - e*N) - s*(d*M + e*C),
			-a*(c*G + s*A) - b*(c*E - s*H),
			a*(c*(d*δγδϕδθ + e*E) + s*(d*D - e*(I + δϕδθ²))) + b*(c*(d*F - e*(N + δϕ²δθ)) - s*(d*(M + δγδϕ²) + e*Z)) - c*(d*C - e*M) + s*(d*N + e*F),
			a*(c*(d*E - e*δγδϕδθ) - s*(d*(I + δϕδθ²) + e*D)) - b*(c*(d*(N + δϕ²δθ) + e*F) + s*(d*Z - e*(M + δγδϕ²))) + c*(d*M + e*C) + s*(d*F - e*N),
			-a*(δ³ϕ - δϕ³) + b*δϕδ²ϕ,
			-a*(d*L + e*B) - b*(d*D - e*I),
			-a*(d*B - e*L) + b*(d*I + e*D),
		)
	elseif rot_seq == :XYX
		return DCM(
			-a*δϕδ²ϕ - b*(δ³ϕ - δϕ³),
			a*(d*D - e*I) - b*(d*L + e*B),
			-a*(d*I + e*D) - b*(d*B - e*L),
			a*(c*E - s*H) - b*(c*G + s*A),
			-a*(c*(d*F - e*(N + δϕ²δθ)) - s*(d*(M + δγδϕ²) + e*Z)) + b*(c*(d*δγδϕδθ + e*E) + s*(d*D - e*(I + δϕδθ²))) - c*(d*C - e*M) + s*(d*N + e*F),
			a*(c*(d*(N + δϕ²δθ) + e*F) + s*(d*Z - e*(M + δγδϕ²))) + b*(c*(d*E - e*δγδϕδθ) - s*(d*(I + δϕδθ²) + e*D)) + c*(d*M + e*C) + s*(d*F - e*N),
			a*(c*H + s*E) + b*(c*A - s*G),
			-a*(c*(d*(M + δγδϕ²) + e*Z) + s*(d*F - e*(N + δϕ²δθ))) - b*(c*(d*D - e*(I + δϕδθ²)) - s*(d*δγδϕδθ + e*E)) - c*(d*N + e*F) - s*(d*C - e*M),
			-a*(c*(d*Z - e*(M + δγδϕ²)) - s*(d*(N + δϕ²δθ) + e*F)) + b*(c*(d*(I + δϕδθ²) + e*D) + s*(d*E - e*δγδϕδθ)) - c*(d*F - e*N) + s*(d*M + e*C),
		)
	elseif rot_seq == :XYZ
		return DCM(
			-a*(d*B - e*L) + b*(d*I + e*D),
			a*(d*L + e*B) + b*(d*D - e*I),
			a*(δ³ϕ - δϕ³) - b*δϕδ²ϕ,
			a*(c*(d*E - e*δγδϕδθ) - s*(d*(I + δϕδθ²) + e*D)) - b*(c*(d*(N + δϕ²δθ) + e*F) + s*(d*Z - e*(M + δγδϕ²))) - c*(d*M + e*C) - s*(d*F - e*N),
			-a*(c*(d*δγδϕδθ + e*E) + s*(d*D - e*(I + δϕδθ²))) - b*(c*(d*F - e*(N + δϕ²δθ)) - s*(d*(M + δγδϕ²) + e*Z)) - c*(d*C - e*M) + s*(d*N + e*F),
			a*(c*G + s*A) + b*(c*E - s*H),
			a*(c*(d*(I + δϕδθ²) + e*D) + s*(d*E - e*δγδϕδθ)) + b*(c*(d*Z - e*(M + δγδϕ²)) - s*(d*(N + δϕ²δθ) + e*F)) + c*(d*F - e*N) - s*(d*M + e*C),
			a*(c*(d*D - e*(I + δϕδθ²)) - s*(d*δγδϕδθ + e*E)) - b*(c*(d*(M + δγδϕ²) + e*Z) + s*(d*F - e*(N + δϕ²δθ))) - c*(d*N + e*F) - s*(d*C - e*M),
			-a*(c*A - s*G) + b*(c*H + s*E),
		)
	elseif rot_seq == :XZX
		return DCM(
			-a*δϕδ²ϕ - b*(δ³ϕ - δϕ³),
			a*(d*I + e*D) + b*(d*B - e*L),
			a*(d*D - e*I) - b*(d*L + e*B),
			-a*(c*H + s*E) - b*(c*A - s*G),
			-a*(c*(d*Z - e*(M + δγδϕ²)) - s*(d*(N + δϕ²δθ) + e*F)) + b*(c*(d*(I + δϕδθ²) + e*D) + s*(d*E - e*δγδϕδθ)) - c*(d*F - e*N) + s*(d*M + e*C),
			a*(c*(d*(M + δγδϕ²) + e*Z) + s*(d*F - e*(N + δϕ²δθ))) + b*(c*(d*D - e*(I + δϕδθ²)) - s*(d*δγδϕδθ + e*E)) + c*(d*N + e*F) + s*(d*C - e*M),
			a*(c*E - s*H) - b*(c*G + s*A),
			-a*(c*(d*(N + δϕ²δθ) + e*F) + s*(d*Z - e*(M + δγδϕ²))) - b*(c*(d*E - e*δγδϕδθ) - s*(d*(I + δϕδθ²) + e*D)) - c*(d*M + e*C) - s*(d*F - e*N),
			-a*(c*(d*F - e*(N + δϕ²δθ)) - s*(d*(M + δγδϕ²) + e*Z)) + b*(c*(d*δγδϕδθ + e*E) + s*(d*D - e*(I + δϕδθ²))) - c*(d*C - e*M) + s*(d*N + e*F),
		)
	elseif rot_seq == :XZY
		return DCM(
			-a*(d*B - e*L) + b*(d*I + e*D),
			-a*(δ³ϕ - δϕ³) + b*δϕδ²ϕ,
			-a*(d*L + e*B) - b*(d*D - e*I),
			-a*(c*(d*(I + δϕδθ²) + e*D) + s*(d*E - e*δγδϕδθ)) - b*(c*(d*Z - e*(M + δγδϕ²)) - s*(d*(N + δϕ²δθ) + e*F)) + c*(d*F - e*N) - s*(d*M + e*C),
			-a*(c*A - s*G) + b*(c*H + s*E),
			a*(c*(d*D - e*(I + δϕδθ²)) - s*(d*δγδϕδθ + e*E)) - b*(c*(d*(M + δγδϕ²) + e*Z) + s*(d*F - e*(N + δϕ²δθ))) + c*(d*N + e*F) + s*(d*C - e*M),
			a*(c*(d*E - e*δγδϕδθ) - s*(d*(I + δϕδθ²) + e*D)) - b*(c*(d*(N + δϕ²δθ) + e*F) + s*(d*Z - e*(M + δγδϕ²))) + c*(d*M + e*C) + s*(d*F - e*N),
			-a*(c*G + s*A) - b*(c*E - s*H),
			a*(c*(d*δγδϕδθ + e*E) + s*(d*D - e*(I + δϕδθ²))) + b*(c*(d*F - e*(N + δϕ²δθ)) - s*(d*(M + δγδϕ²) + e*Z)) - c*(d*C - e*M) + s*(d*N + e*F),
		)
	elseif rot_seq == :YXY
		return DCM(
			-a*(c*(d*F - e*(N + δϕ²δθ)) - s*(d*(M + δγδϕ²) + e*Z)) + b*(c*(d*δγδϕδθ + e*E) + s*(d*D - e*(I + δϕδθ²))) - c*(d*C - e*M) + s*(d*N + e*F),
			a*(c*E - s*H) - b*(c*G + s*A),
			-a*(c*(d*(N + δϕ²δθ) + e*F) + s*(d*Z - e*(M + δγδϕ²))) - b*(c*(d*E - e*δγδϕδθ) - s*(d*(I + δϕδθ²) + e*D)) - c*(d*M + e*C) - s*(d*F - e*N),
			a*(d*D - e*I) - b*(d*L + e*B),
			-a*δϕδ²ϕ - b*(δ³ϕ - δϕ³),
			a*(d*I + e*D) + b*(d*B - e*L),
			a*(c*(d*(M + δγδϕ²) + e*Z) + s*(d*F - e*(N + δϕ²δθ))) + b*(c*(d*D - e*(I + δϕδθ²)) - s*(d*δγδϕδθ + e*E)) + c*(d*N + e*F) + s*(d*C - e*M),
			-a*(c*H + s*E) - b*(c*A - s*G),
			-a*(c*(d*Z - e*(M + δγδϕ²)) - s*(d*(N + δϕ²δθ) + e*F)) + b*(c*(d*(I + δϕδθ²) + e*D) + s*(d*E - e*δγδϕδθ)) - c*(d*F - e*N) + s*(d*M + e*C),
		)
	elseif rot_seq == :YXZ
		return DCM(
			a*(c*(d*δγδϕδθ + e*E) + s*(d*D - e*(I + δϕδθ²))) + b*(c*(d*F - e*(N + δϕ²δθ)) - s*(d*(M + δγδϕ²) + e*Z)) - c*(d*C - e*M) + s*(d*N + e*F),
			a*(c*(d*E - e*δγδϕδθ) - s*(d*(I + δϕδθ²) + e*D)) - b*(c*(d*(N + δϕ²δθ) + e*F) + s*(d*Z - e*(M + δγδϕ²))) + c*(d*M + e*C) + s*(d*F - e*N),
			-a*(c*G + s*A) - b*(c*E - s*H),
			-a*(d*L + e*B) - b*(d*D - e*I),
			-a*(d*B - e*L) + b*(d*I + e*D),
			-a*(δ³ϕ - δϕ³) + b*δϕδ²ϕ,
			a*(c*(d*D - e*(I + δϕδθ²)) - s*(d*δγδϕδθ + e*E)) - b*(c*(d*(M + δγδϕ²) + e*Z) + s*(d*F - e*(N + δϕ²δθ))) + c*(d*N + e*F) + s*(d*C - e*M),
			-a*(c*(d*(I + δϕδθ²) + e*D) + s*(d*E - e*δγδϕδθ)) - b*(c*(d*Z - e*(M + δγδϕ²)) - s*(d*(N + δϕ²δθ) + e*F)) + c*(d*F - e*N) - s*(d*M + e*C),
			-a*(c*A - s*G) + b*(c*H + s*E),
		)
	elseif rot_seq == :YZX
		return DCM(
			-a*(c*A - s*G) + b*(c*H + s*E),
			a*(c*(d*(I + δϕδθ²) + e*D) + s*(d*E - e*δγδϕδθ)) + b*(c*(d*Z - e*(M + δγδϕ²)) - s*(d*(N + δϕ²δθ) + e*F)) + c*(d*F - e*N) - s*(d*M + e*C),
			a*(c*(d*D - e*(I + δϕδθ²)) - s*(d*δγδϕδθ + e*E)) - b*(c*(d*(M + δγδϕ²) + e*Z) + s*(d*F - e*(N + δϕ²δθ))) - c*(d*N + e*F) - s*(d*C - e*M),
			a*(δ³ϕ - δϕ³) - b*δϕδ²ϕ,
			-a*(d*B - e*L) + b*(d*I + e*D),
			a*(d*L + e*B) + b*(d*D - e*I),
			a*(c*G + s*A) + b*(c*E - s*H),
			a*(c*(d*E - e*δγδϕδθ) - s*(d*(I + δϕδθ²) + e*D)) - b*(c*(d*(N + δϕ²δθ) + e*F) + s*(d*Z - e*(M + δγδϕ²))) - c*(d*M + e*C) - s*(d*F - e*N),
			-a*(c*(d*δγδϕδθ + e*E) + s*(d*D - e*(I + δϕδθ²))) - b*(c*(d*F - e*(N + δϕ²δθ)) - s*(d*(M + δγδϕ²) + e*Z)) - c*(d*C - e*M) + s*(d*N + e*F),
		)
	elseif rot_seq == :YZY
		return DCM(
			-a*(c*(d*Z - e*(M + δγδϕ²)) - s*(d*(N + δϕ²δθ) + e*F)) + b*(c*(d*(I + δϕδθ²) + e*D) + s*(d*E - e*δγδϕδθ)) - c*(d*F - e*N) + s*(d*M + e*C),
			a*(c*H + s*E) + b*(c*A - s*G),
			-a*(c*(d*(M + δγδϕ²) + e*Z) + s*(d*F - e*(N + δϕ²δθ))) - b*(c*(d*D - e*(I + δϕδθ²)) - s*(d*δγδϕδθ + e*E)) - c*(d*N + e*F) - s*(d*C - e*M),
			-a*(d*I + e*D) - b*(d*B - e*L),
			-a*δϕδ²ϕ - b*(δ³ϕ - δϕ³),
			a*(d*D - e*I) - b*(d*L + e*B),
			a*(c*(d*(N + δϕ²δθ) + e*F) + s*(d*Z - e*(M + δγδϕ²))) + b*(c*(d*E - e*δγδϕδθ) - s*(d*(I + δϕδθ²) + e*D)) + c*(d*M + e*C) + s*(d*F - e*N),
			a*(c*E - s*H) - b*(c*G + s*A),
			-a*(c*(d*F - e*(N + δϕ²δθ)) - s*(d*(M + δγδϕ²) + e*Z)) + b*(c*(d*δγδϕδθ + e*E) + s*(d*D - e*(I + δϕδθ²))) - c*(d*C - e*M) + s*(d*N + e*F),
		)
	elseif rot_seq == :ZXY
		return DCM(
			-a*(c*(d*δγδϕδθ + e*E) + s*(d*D - e*(I + δϕδθ²))) - b*(c*(d*F - e*(N + δϕ²δθ)) - s*(d*(M + δγδϕ²) + e*Z)) - c*(d*C - e*M) + s*(d*N + e*F),
			a*(c*G + s*A) + b*(c*E - s*H),
			a*(c*(d*E - e*δγδϕδθ) - s*(d*(I + δϕδθ²) + e*D)) - b*(c*(d*(N + δϕ²δθ) + e*F) + s*(d*Z - e*(M + δγδϕ²))) - c*(d*M + e*C) - s*(d*F - e*N),
			a*(c*(d*D - e*(I + δϕδθ²)) - s*(d*δγδϕδθ + e*E)) - b*(c*(d*(M + δγδϕ²) + e*Z) + s*(d*F - e*(N + δϕ²δθ))) - c*(d*N + e*F) - s*(d*C - e*M),
			-a*(c*A - s*G) + b*(c*H + s*E),
			a*(c*(d*(I + δϕδθ²) + e*D) + s*(d*E - e*δγδϕδθ)) + b*(c*(d*Z - e*(M + δγδϕ²)) - s*(d*(N + δϕ²δθ) + e*F)) + c*(d*F - e*N) - s*(d*M + e*C),
			a*(d*L + e*B) + b*(d*D - e*I),
			a*(δ³ϕ - δϕ³) - b*δϕδ²ϕ,
			-a*(d*B - e*L) + b*(d*I + e*D),
		)
	elseif rot_seq == :ZXZ
		return DCM(
			-a*(c*(d*F - e*(N + δϕ²δθ)) - s*(d*(M + δγδϕ²) + e*Z)) + b*(c*(d*δγδϕδθ + e*E) + s*(d*D - e*(I + δϕδθ²))) - c*(d*C - e*M) + s*(d*N + e*F),
			a*(c*(d*(N + δϕ²δθ) + e*F) + s*(d*Z - e*(M + δγδϕ²))) + b*(c*(d*E - e*δγδϕδθ) - s*(d*(I + δϕδθ²) + e*D)) + c*(d*M + e*C) + s*(d*F - e*N),
			a*(c*E - s*H) - b*(c*G + s*A),
			-a*(c*(d*(M + δγδϕ²) + e*Z) + s*(d*F - e*(N + δϕ²δθ))) - b*(c*(d*D - e*(I + δϕδθ²)) - s*(d*δγδϕδθ + e*E)) - c*(d*N + e*F) - s*(d*C - e*M),
			-a*(c*(d*Z - e*(M + δγδϕ²)) - s*(d*(N + δϕ²δθ) + e*F)) + b*(c*(d*(I + δϕδθ²) + e*D) + s*(d*E - e*δγδϕδθ)) - c*(d*F - e*N) + s*(d*M + e*C),
			a*(c*H + s*E) + b*(c*A - s*G),
			a*(d*D - e*I) - b*(d*L + e*B),
			-a*(d*I + e*D) - b*(d*B - e*L),
			-a*δϕδ²ϕ - b*(δ³ϕ - δϕ³),
		)
	elseif rot_seq == :ZYZ
		return DCM(
			-a*(c*(d*Z - e*(M + δγδϕ²)) - s*(d*(N + δϕ²δθ) + e*F)) + b*(c*(d*(I + δϕδθ²) + e*D) + s*(d*E - e*δγδϕδθ)) - c*(d*F - e*N) + s*(d*M + e*C),
			a*(c*(d*(M + δγδϕ²) + e*Z) + s*(d*F - e*(N + δϕ²δθ))) + b*(c*(d*D - e*(I + δϕδθ²)) - s*(d*δγδϕδθ + e*E)) + c*(d*N + e*F) + s*(d*C - e*M),
			-a*(c*H + s*E) - b*(c*A - s*G),
			-a*(c*(d*(N + δϕ²δθ) + e*F) + s*(d*Z - e*(M + δγδϕ²))) - b*(c*(d*E - e*δγδϕδθ) - s*(d*(I + δϕδθ²) + e*D)) - c*(d*M + e*C) - s*(d*F - e*N),
			-a*(c*(d*F - e*(N + δϕ²δθ)) - s*(d*(M + δγδϕ²) + e*Z)) + b*(c*(d*δγδϕδθ + e*E) + s*(d*D - e*(I + δϕδθ²))) - c*(d*C - e*M) + s*(d*N + e*F),
			a*(c*E - s*H) - b*(c*G + s*A),
			a*(d*I + e*D) + b*(d*B - e*L),
			a*(d*D - e*I) - b*(d*L + e*B),
			-a*δϕδ²ϕ - b*(δ³ϕ - δϕ³),
		)
	else
		throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
	end

end
