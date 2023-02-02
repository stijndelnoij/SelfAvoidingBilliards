export collision, next_collision

#####################################################################################
# Accuracy & Convenience Functions
#####################################################################################
"""
Approximate arccos(1 - x) for x very close to 0.
"""
@inline (acos1mx(x::T)::T) where {T} = sqrt(2x) + sqrt(x)^3/sixsqrt
const sixsqrt = 6sqrt(2)

@inline cross2D(a, b) = a[1]*b[2] - a[2]*b[1]

@inline accuracy(::Type{T}) where {T} = sqrt(eps(T))
@inline accuracy(::Type{BigFloat}) = BigFloat(1e-32)
@inline nocollision(::Type{T}) where {T} = (T(Inf), SV{T}(0.0, 0.0))
@inline noresolution(::Type{T}) where {T} = (T(Inf), SV{T}(0.0, 0.0))

#######################################################################################
## Particle
#######################################################################################
"""
    collision(p::AbstractParticle, o::Obstacle) → t, cp
Find the collision (if any) between given particle and obstacle.
Return the time until collision and the estimated collision point `cp`.

Return `Inf, SV(0, 0)` if the collision is not possible *or* if the
collision happens backwards in time.

**It is the duty of `collision` to avoid incorrect collisions when the particle is
on top of the obstacle (or very close).**
"""
function collision(p::Particle{T}, w::Wall{T}) where {T}
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    if denom ≥ 0.0
        return nocollision(T)
    else
        t = dot(w.sp - p.pos, n)/denom
        return t, p.pos + t*p.vel
    end
end

function collision(p::Particle{T}, o::Union{FiniteWall{T}, Tail{T}}) where {T}
    vdelta = [[o.ep[1]-o.sp[1], o.ep[2]-o.sp[2]] [-p.vel[1], -p.vel[2]]] #Matrix with segment- and velocity vector
    vdeltadet = abs(det(vdelta)) #A meassure for how close the particle and segment are
    if vdeltadet > accuracy(T) #If the segment and particle are to close we have encountered a tunneling event.
                                    #The collision is false in case of a tunneling event. So the particle 
                                    #ignores the collision and moves like it would have without tunneling.
        tvector = inv(vdelta)*(p.pos-o.sp) #Vector with collision time inside segment and col-time for the particle.
        t = tvector[2]
        
        if 0 < tvector[1] < 1 && tvector[2] > accuracy(T) # Collision must be on-line and above minimal travel time.
            return t, p.pos + t * p.vel
        else
            return nocollision(T)
        end
        
    else 
        return noresolution(T)
    end
end
function collision(p::Particle{T}, d::Circular{T}) where {T}

    dotp = dot(p.vel, normalvec(d, p.pos))
    dotp ≥ 0.0 && return nocollision(T)

    dc = p.pos - d.c
    B = dot(p.vel, dc)           #pointing towards circle center: B < 0
    C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
    Δ = B*B - C

    Δ ≤ 0.0 && return nocollision(T)
    sqrtD = sqrt(Δ)

    # Closest point:
    t = -B - sqrtD
    return t, p.pos + t*p.vel
end

function collision(p::Particle{T}, d::Antidot{T}) where {T}

    dotp = dot(p.vel, normalvec(d, p.pos))
    if d.pflag == true
        dotp ≥ 0 && return nocollision(T)
    end

    dc = p.pos - d.c
    B = dot(p.vel, dc)           #velocity towards circle center: B < 0
    C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
    Δ = B^2 - C
    Δ ≤ 0 && return nocollision(T)
    sqrtD = sqrt(Δ)

    # Closest point (may be in negative time):
    if dotp < 0
        t = -B - sqrtD
    else
        t = -B + sqrtD
    end
    t ≤ 0.0 ? nocollision(T) : (t, p.pos + t*p.vel)
end

function collision(p::Particle{T}, d::Semicircle{T}) where {T}

    dc = p.pos - d.c
    B = dot(p.vel, dc)         #velocity towards circle center: B > 0
    C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
    Δ = B^2 - C

    Δ ≤ 0 && return nocollision(T)
    sqrtD = sqrt(Δ)

    nn = dot(dc, d.facedir)
    if nn ≥ 0 # I am NOT inside semicircle
        # Return most positive time
        t = -B + sqrtD
    else # I am inside semicircle:
        t = -B - sqrtD
        # these lines make sure that the code works for ANY starting position:
        if t ≤ 0 || distance(p, d) ≤ accuracy(T)
            t = -B + sqrtD
        end
    end
    # This check is necessary to not collide with the non-existing side
    newpos = p.pos + p.vel * t
    if dot(newpos - d.c, d.facedir) ≥ 0 # collision point on BAD HALF;
        return nocollision(T)
    end
    # If collision time is negative, return Inf:
    t ≤ 0.0 ? nocollision(T) : (t, p.pos + t*p.vel)
end


function collision(p::Particle{T}, e::Ellipse{T}) where {T}
    # First check if particle is "looking at" eclipse if it is outside
    if e.pflag
        # These lines may be "not accurate enough" but so far all is good
        dotp = dot(p.vel, normalvec(e, p.pos))
        dotp ≥ 0.0 && return nocollision(T)
    end

    # http://www.ambrsoft.com/TrigoCalc/Circles2/Ellipse/EllipseLine.htm
    a = e.a; b = e.b
    # Translate particle with ellipse center (so that ellipse lies on [0, 0])
    pc = p.pos - e.c
    # Find μ, ψ for line equation y = μx + ψ describing particle
    μ = p.vel[2]/p.vel[1]
    ψ = pc[2] - μ*pc[1]

    # Determinant and intersection points follow from the link
    denomin = a*a*μ*μ + b*b
    Δ² = denomin - ψ*ψ
    Δ² ≤ 0 && return nocollision(T)
    Δ = sqrt(Δ²); f1 = -a*a*μ*ψ; f2 = b*b*ψ # just factors
    I1 = SV(f1 + a*b*Δ, f2 + a*b*μ*Δ)/denomin
    I2 = SV(f1 - a*b*Δ, f2 - a*b*μ*Δ)/denomin

    d1 = norm(pc - I1); d2 = norm(pc - I2)
    if e.pflag
        return d1 < d2 ? (d1, I1 + e.c) : (d2, I2 + e.c)
    else # inside the ellipse: one collision is _always_ valid
        if d1 < d2
            dmin, Imin = d1, I1
            dmax, Imax = d2, I2
        else
            dmin, Imin = d2, I2
            dmax, Imax = d1, I1
        end

        if dmin < accuracy(T) # special case for being very close to ellipse
            dotp = dot(p.vel, normalvec(e, Imin))
            dotp ≥ 0 && return (dmax, Imax + e.c)
        end
         # check which of the two points is ahead or behind the obstacle
        z1 = dot(pc - Imax, p.vel)
        return z1 < 0 ? (dmax, Imax + e.c) : (dmin, Imin + e.c)
    end
end



#######################################################################################
## next_collision
#######################################################################################
"""
    next_collision(p::AbstractParticle, bd::Billiard) -> i, tmin, cp
Compute the [`collision`](@ref) across all obstacles in `bd` and find the minimum
one. Return the index of colliding obstacle, the time and the collision point.
"""

function next_collision(p, bd::Billiard{T}) where {T}
    ind = 0
    tmin = T(Inf)
    cp = SV{T}(0.0, 0.0)
    for j in 1:length(bd)
        x = bd[j]
        tcol, pcol = collision(p, x)
        if tcol == -T(Inf)
            tmin = tcol
            ind = 0
            cp = pcol
            break
        end
        # Set minimum time:
        if tcol < tmin
            tmin = tcol
            ind = j
            cp = pcol
        end
    end            
    return ind, tmin, cp
end
