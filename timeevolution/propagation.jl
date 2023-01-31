export bounce!, ispinned, propagate!

@inline increment_counter(::Int, t_to_write) = 1
@inline increment_counter(::T, t_to_write) where {T<:AbstractFloat} = t_to_write

#####################################################################################
# Bounce
#####################################################################################
"""
    bounce!(p::AbstractParticle, bd::Billiard) → i, t, pos, vel
"Bounce" the particle (advance for one collision) in the billiard.
Takes care of finite-precision issues.

Return:
* index of the obstacle that the particle just collided with
* the time from the previous collision until the current collision `t`
* position and velocity of the particle at the current collision (*after* the
  collision has been resolved!). The position is given in the unit cell of
  periodic billiards. Do `pos += p.current_cell` for the position in real space.

```julia
bounce!(p, bd, raysplit) → i, t, pos, vel
```
Ray-splitting version of `bounce!`.
"""
@inline function bounce!(p::AbstractParticle{T}, bd::Billiard{T}, collisions::DataFrame) where {T}
    i::Int, tmin::T, cp::SV{T} = next_collision(p, bd)
    α = 0.0
    if tmin != T(Inf)
        o = bd[i]
        vel = cp - p.pos
        w = cp - p.pos
        normal = [-w[2], w[1]]
        dotpr = abs(dot(vel, o.normal))/(norm(vel)*norm(o.normal))
        dotpr <= 1. ? α = acos(dotpr) : α = 0.0
    end
    if α != 0.0 && tmin != T(Inf)
        new_wall_id = last(bd).id + 1
        new_wall = Tail(p.pos, cp, normal, new_wall_id, new_wall_id + 1)
        relocate!(p, o, tmin, cp)
        resolvecollision!(p, o)
        col = add_collision(o, new_wall, vel)
        push!(collisions, col)
        push!(bd.obstacles, new_wall)
    end
    return i, tmin, p.pos, p.vel
end


#####################################################################################
# Relocate
#####################################################################################
"""
    relocate!(p::AbstractParticle, o::Obstacle, t, cp)
Propagate the particle to `cp` and propagate velocities for time `t`.
Check if it is on the correct side of the obstacle. If not,
change the particle position by [`distance`](@ref) along the [`normalvec`](@ref)
of the obstacle.
"""
function relocate!(p::AbstractParticle{T}, o::Obstacle{T}, tmin::T, cp::SV{T}) where {T}
    propagate!(p, cp, tmin) # propagate to collision point
    okay = _correct_pos!(p, o)
    return okay
end
function _correct_pos!(p, o)
    d = distance(p.pos, o)
    okay = _okay(d, o)
    if !okay
        n = normalvec(o, p.pos)
        p.pos -= d * n
        # due to finite precision this does not always give positive distance
        # but collision takes care of it, as it does not care about collisions
        # very close to the point.
    end
    return okay
end

_okay(d, o::Obstacle) = d ≥ 0.0
_okay(d, o::PeriodicWall) = d ≤ 0.0

"""
    propagate!(p::AbstractParticle, t)
Propagate the particle `p` for given time `t`, changing appropriately the the
`p.pos` and `p.vel` fields.

    propagate!(p, position, t)
Do the same, but take advantage of the already calculated `position` that the
particle should end up at.
"""
@inline function propagate!(p::Particle{T}, t::Real) where {T}
    p.pos += SV{T}(p.vel[1]*t, p.vel[2]*t)
end
@inline propagate!(p::Particle, newpos::SV, t::Real) = (p.pos = newpos)



#####################################################################################
# Resolve Collisions
#####################################################################################
"""
    specular!(p::AbstractParticle, o::Obstacle)
Perform specular reflection based on the normal vector of the Obstacle.

In the case where the given obstacle is a `RandomObstacle`, the specular reflection
randomizes the velocity instead (within -π/2+ε to π/2-ε of the normal vector).
"""
@inline function specular!(p::AbstractParticle{T}, o::Obstacle{T})::Nothing where {T}
    n = normalvec(o, p.pos)
    p.vel = p.vel - 2*dot(n, p.vel)*n
    return nothing
end

@inline function specular!(p::AbstractParticle{T}, r::RandomDisk{T})::Nothing where {T}
    n = normalvec(r, p.pos)
    φ = atan(n[2], n[1]) + 0.95(π*rand() - π/2) #this cannot be exactly π/2
    p.vel = SVector{2,T}(cos(φ), sin(φ))
    return nothing
end

@inline function specular!(p::AbstractParticle{T}, r::RandomWall{T})::Nothing where {T}
    n = normalvec(r, p.pos)
    φ = atan(n[2], n[1]) + 0.95(π*rand() - π/2) #this cannot be exactly π/2
    p.vel = SVector{2,T}(cossin(φ))
    return nothing
end

"""
    periodicity!(p::AbstractParticle, w::PeriodicWall)
Perform periodicity conditions of `w` on `p`.
"""
@inline function periodicity!(p::AbstractParticle, w::PeriodicWall)::Nothing
    p.pos += w.normal
    p.current_cell -= w.normal
    return nothing
end

"""
    resolvecollision!(p::AbstractParticle, o::Obstacle)
Resolve the collision between particle `p` and obstacle `o`, depending on the
type of `o` (do `specular!` or `periodicity!`).

    resolvecollision!(p, o, T::Function, θ::Function, new_ω::Function)
This is the ray-splitting implementation. The three functions given are drawn from
the ray-splitting dictionary that is passed directly to `evolve!()`. For a calculated
incidence angle φ, if T(φ) > rand(), ray-splitting occurs.
"""
@inline resolvecollision!(p::Particle, o::Obstacle) = specular!(p, o)
@inline resolvecollision!(p::Particle, o::PeriodicWall) = periodicity!(p, o)
