export AbstractParticle, Particle
####################################################
## Particle
####################################################
"""
    AbstractParticle
Particle supertype.
"""
abstract type AbstractParticle{T<:AbstractFloat} end
eltype(p::AbstractParticle{T}) where {T} = T


"""
```julia
Particle(ic::Vector{T}) #where ic = [x0, y0, φ0]
Particle(x0, y0, φ0)
Particle(pos::SVector, vel::SVector)
```
Create a particle with initial conditions `x0, y0, φ0`. It propagates as
a straight line.

The field `current_cell` shows at which cell of a periodic billiard is the particle
currently located.
"""
mutable struct Particle{T<:AbstractFloat} <: AbstractParticle{T}
    pos::SVector{2,T}
    vel::SVector{2,T}
    current_cell::SVector{2,T}
    function Particle(
        pos::SVector{2,T}, vel::SVector{2,T}, cc::SVector{2,T}) where{T<:AbstractFloat}
        new{T}(pos, normalize(vel), cc)
    end
end

Base.copy(p::Particle) =
Particle(p.pos, p.vel, p.current_cell)

function Particle(ic::AbstractVector{S}) where {S<:Real}
    T = S<:Integer ? Float64 : S
    φ0 = ic[3]
    pos = SVector{2,T}(ic[1:2]); vel = SVector{2,T}(cossin(φ0)...)
    return Particle(pos, vel, SVector{2,T}(0,0))
end
Particle(x::Real, y::Real, φ::Real) = Particle(collect(promote(x,y,φ)))
Particle() = Particle(rand(), rand(), rand()*2π)
function Particle(pos::SV{T}, vel::SV{T}) where {T}
    S = T<:Integer ? Float64 : T
    return Particle(pos, vel, SVector{2,S}(0.0, 0.0))
end
show(io::IO, p::Particle{T}) where {T} =
print(io, "Particle{$T}\n",
"position: $(p.pos+p.current_cell)\nvelocity: $(p.vel)")

