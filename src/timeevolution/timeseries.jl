export evolve!, evolve, timeseries!, timeseries

"""
    evolve!([p::AbstractParticle,] bd::Billiard, t, df::DataFrame)
Evolve the given particle `p` inside the billiard `bd`. If `t` is of type
`AbstractFloat`, evolve for as much time as `t`. If however `t` is of type `Int`,
evolve for as many collisions as `t`.
Return the states of the particle between collisions.

This function mutates the particle, use `evolve` otherwise. If a particle is
not given, a random one is picked through [`randominside`](@ref).

### Return

* `ct::Vector{T}` : Collision times.
* `poss::Vector{SVector{2,T}}` : Positions at the collisions.
* `vels::Vector{SVector{2,T}})` : Velocities exactly after the collisions.

The time `ct[i+1]` is the time necessary to reach state `poss[i+1], vels[i+1]` starting
from the state `poss[i], vels[i]`. That is why `ct[1]` is always 0 since
`poss[1], vels[1]` are the initial conditions.

Notice that at any point, the velocity vector `vels[i]` is the one obdained *after*
the specular reflection of the `i-1`th collision.
"""
function evolve!(p::AbstractParticle{T}, bd::Billiard{T}, t, df::DataFrame;
    warning = false) where {T<:AbstractFloat}

    if t ≤ 0
        throw(ArgumentError("`evolve!()` cannot evolve backwards in time."))
    end

    rt = T[0.0]; rpos = [p.pos]; rvel = [p.vel]
    count = zero(t)
    flight = zero(T)
    if typeof(t) == Int
        for zzz in (rt, rpos, rvel)
            sizehint!(zzz, t)
        end
    end

    while count < t

        i, tmin, pos, vel = bounce!(p, bd, df)
        flight += tmin

        if isperiodic(i, bd)
            continue
        else
            push!(rpos, pos + p.current_cell)
            push!(rvel, vel)
            push!(rt, flight)
            # set counter
            count += increment_counter(t, flight)
            flight = zero(T)
        end
    end#time, or collision number, loop

    # Return stuff
    if isray
        return (rt, rpos, rvel, omegas)
    else
        return (rt, rpos, rvel)
    end
end

"""
    evolve(p, args...)
Same as [`evolve!`](@ref) but copies the particle instead.
"""
evolve(p::AbstractParticle, args...) = evolve!(copy(p), args...)
evolve(bd::Billiard, args...; kwargs...) =
    evolve!(randominside(bd), bd, args...; kwargs...)

#####################################################################################
# Timeseries
#####################################################################################
"""
    timeseries!([p::AbstractParticle,] bd::Billiard, t; dt, warning)
Evolve the given particle `p` inside the billiard `bd` for the condition `t`
and return the x, y, vx, vy timeseries and the time vector.
If `t` is of type `AbstractFloat`, then evolve for as much time as `t`.
If however `t` is of type `Int`, evolve for as many collisions as `t`.
Otherwise, `t` can be any function, that takes as an input `t(n, τ, i, p)`
and returns `true` when the evolution should terminate. Here `n` is the amount
of obstacles collided with so far, `τ` the amount time evolved so far, `i`
the obstacle just collided with and `p` the particle (so you can access e.g.
`p.pos`).

This function mutates the particle, use `timeseries` otherwise. If a particle is
not given, a random one is picked through [`randominside`](@ref).

The keyword argument `dt` is the time step used for interpolating the time
series in between collisions. `dt` is capped by the collision time, as
the interpolation _always_ stops at collisions.
For straight propagation `dt = Inf`.

Internally uses [`DynamicalBilliards.extrapolate`](@ref).
"""
function timeseries!(p::AbstractParticle{T}, bd::Billiard{T}, f, df::DataFrame;
                     dt = typeof(p) <: Particle ? T(Inf) : T(0.01),
                     warning::Bool = true) where {T}

    ts = [zero(T)]
    x  = [p.pos[1]]; y  = [p.pos[2]]
    vx = [p.vel[1]]; vy = [p.vel[2]]
    n, i, t, flight = 0, 0, zero(T), zero(T)
    prevpos = p.pos + p.current_cell
    prevvel = p.vel

    @inbounds while check_condition(f, n, t, i, p) # count < t
        i, ct = bounce!(p, bd, df)
        flight += ct
        if isperiodic(i, bd) # do nothing at periodic obstacles
            continue
        else
            if flight ≤ dt # push collision point only
                push!(ts, flight + t)
                push!(x, p.pos[1] + p.current_cell[1])
                push!(y, p.pos[2] + p.current_cell[2])
                push!(vx, p.vel[1]); push!(vy, p.vel[2])
            else
                nx, ny, nvx, nvy, nts = extrapolate(p, prevpos, prevvel, flight, dt)
                append!(ts, nts[2:end] .+ t)
                append!(x, nx[2:end]); append!(vx, nvx[2:end])
                append!(y, ny[2:end]); append!(vy, nvy[2:end])
            end
            prevpos = p.pos + p.current_cell
            prevvel = p.vel
            t += flight; n += 1
            flight = zero(T)
        end
    end
    return x, y, vx, vy, ts
end

check_condition(f::Int, n, t, i, p) = f > n
check_condition(f::AbstractFloat, n, t, i, p) = f > t
check_condition(f, n, t, i, p) = !f(n, t, i, p)

"""
    extrapolate(particle, prevpos, prevvel, ct, dt) → x, y, vx, vy, t

Create the timeseries that connect a `particle`'s previous position and velocity
`prevpos, prevvel` with the `particle`'s current position and velocity,
provided that the collision time between previous and current state is `ct`.

`dt` is the sampling time, the angular velocity that governed the free flight.

Here is how this function is used (for example)
```julia
prevpos, prevvel = p.pos + p.current_cell, p.vel
for _ in 1:5
    i, ct, pos, vel = bounce!(p, bd, df)
    x, y, x, vy, t = extrapolate(p, prevpos, prevvel, ct, 0.1)
    # append x, y, ... to other vectors or do whatever with them
    prevpos, prevvel = p.pos + p.current_cell, p.vel
end
```
"""

function extrapolate(p::Particle{T}, prevpos::SV{T}, prevvel::SV{T}, ct::T,
                     dt::T) where {T <: AbstractFloat}

    tvec = collect(0:dt:ct)
    x = Vector{T}(undef, length(tvec))
    y = Vector{T}(undef, length(tvec))
    vx = Vector{T}(undef, length(tvec))
    vy = Vector{T}(undef, length(tvec))

    @inbounds for (i,t) ∈ enumerate(tvec)
        x[i] = prevpos[1] + t*prevvel[1]
        y[i] = prevpos[2] + t*prevvel[2]
        vx[i], vy[i] = prevvel
    end

    # finish with ct
    @inbounds if tvec[end] != ct
        push!(tvec, ct)
        push!(x, p.pos[1] + p.current_cell[1])
        push!(y, p.pos[2] + p.current_cell[2])
        push!(vx, p.vel[1]); push!(vy, p.vel[2])
    end

    return x, y, vx, vy, tvec
end

# non-mutating version
"""
    timeseries(p, args...; kwargs...)
Non-mutating version of [`timeseries!`](@ref)
"""
@inline timeseries(p::AbstractParticle, args...; kwargs...) =
    timeseries!(copy(p), args...; kwargs...)
