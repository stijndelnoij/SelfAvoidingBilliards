export Obstacle, Disk, Antidot, RandomDisk, Wall, Circular,
InfiniteWall, PeriodicWall, RandomWall, FiniteWall,
Semicircle, Ellipse, normalvec, prev_obst
export translate

using InteractiveUtils
#######################################################################################
## Circles
#######################################################################################
"""
    Obstacle{<:AbstractFloat}
Obstacle supertype.
"""
abstract type Obstacle{T<:AbstractFloat} end
eltype(o::Obstacle{T}) where {T} = T


"""
    Circular{T<:AbstractFloat} <: Obstacle{T}
Circular obstacle supertype.
"""
abstract type Circular{T<:AbstractFloat} <: Obstacle{T} end

"""
    Disk{T<:AbstractFloat}  <: Circular{T}
Disk-like obstacle with propagation allowed outside of the circle (immutable type).
### Fields:
* `c::SVector{2,T}` : Center.
* `r::T` : Radius.
* `id::Integer` : Should be unique for every obstacle.
"""
struct Disk{T<:AbstractFloat} <: Circular{T}
    c::SVector{2,T}
    r::T
    id::Integer
    next::Integer
end

function Disk(c::AbstractVector{T}, r::Real, id::Integer) where {T<:Real}
    S = T <: Integer ? Float64 : T
    return Disk{S}(SVector{2,S}(c), convert(S, abs(r)), id, id)
end
Disk{T}(args...) where {T} = Disk(args...)

"""
    RandomDisk{T<:AbstractFloat} <: Circular{T}
Disk-like obstacle that randomly (and uniformly) reflects colliding particles.
The propagation is allowed outside of the circle.
### Fields:
* `c::SVector{2,T}` : Center.
* `r::T` : Radius.
* `id::Integer` : Obstacle id.
"""
struct RandomDisk{T<:AbstractFloat} <: Circular{T}
    c::SVector{2,T}
    r::T
    id::Integer
    next::Integer
end
function RandomDisk(
    c::AbstractVector{T}, r::Real, id::Integer) where {T<:Real}
    S = T <: Integer ? Float64 : T
    return RandomDisk{S}(SVector{2,S}(c), convert(S, abs(r)), id, id)
end
RandomDisk{T}(args...) where {T} = RandomDisk(args...)

"""
    Antidot{T<:AbstractFloat} <: Circular{T}
Disk-like obstacle that allows propagation both inside and outside of the disk
(mutable type). Used in ray-splitting billiards.
### Fields:
* `c::SVector{2,T}` : Center.
* `r::T` : Radius.
* `pflag::Bool` : Flag that keeps track of where the particle is currently
  propagating (`pflag` = propagation-flag).
  `true` stands for *outside* the disk, `false` for *inside* the disk.
  Defaults to `true`.
  * `id::Integer` : Obstacle id.
"""
mutable struct Antidot{T<:AbstractFloat} <: Circular{T}
    c::SVector{2,T}
    r::T
    pflag::Bool
    id::Integer
    next::Integer
end
function Antidot(c::AbstractVector{T}, r::Real,
    pflag::Bool = true, id::Integer) where {T<:Real}
    S = T <: Integer ? Float64 : T
    return Antidot{S}(SVector{2,S}(c), convert(S, abs(r)), pflag, id, id)
end
Antidot(c, r, id::Integer) = Antidot(c, r, true, id)
Antidot{T}(args...) where {T} = Antidot(args...)

"""
    Semicircle{T<:AbstractFloat} <: Circular{T}
Obstacle that represents half a circle. Propagation is allowed only inside
the semicircle.
### Fields:
* `c::SVector{2,T}` : Center.
* `r::T` : Radius.
* `facedir::SVector{2,T}` : Direction where the open face of the Semicircle is facing.
* `id::Integer` : Obstacle id.
"""
struct Semicircle{T<:AbstractFloat} <: Circular{T}
    c::SVector{2,T}
    r::T
    facedir::SVector{2,T}
    id::Integer
    next::Integer
end
function Semicircle(
    c::AbstractVector{T}, r::Real, facedir, id::Integer) where {T<:Real}
    S = T <: Integer ? Float64 : T
    return Semicircle{S}(
    SV(c), convert(S, abs(r)), SV(normalize(facedir)), id, id)
end

show(io::IO, w::Circular{T}) where {T} =
print(io, "Circular object $(w.id) {$T}\n", "center: $(w.c)\nradius: $(w.r)")

show(io::IO, w::Semicircle{T}) where {T} =
print(io, "Circular object $(w.id) {$T}\n", "center: $(w.c)\nradius: $(w.r)\nfacedir: $(w.facedir)")



#######################################################################################
## Walls
#######################################################################################
"""
    Wall{T<:AbstractFloat} <: Obstacle{T}
Wall obstacle supertype.

All `Wall` subtypes (except `PeriodicWall`) can be called as `Wall(sp, ep)`,
in which case the normal vector is computed automatically to point to the left
of `v = ep - sp`.
"""
abstract type Wall{T<:AbstractFloat} <: Obstacle{T} end

"""
    InfiniteWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing specular reflection during collision (immutable type).
Faster than [`FiniteWall`](@ref), meant to be used for convex billiards.
### Fields:
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the
  particle *will come from before a collision* (pointing towards the inside of the
  billiard). The size of the vector is irrelevant
  since it is internally normalized.
  * `id::Integer` : Obstacle id.
"""
struct InfiniteWall{T<:AbstractFloat} <: Wall{T}
    sp::SVector{2,T}
    ep::SVector{2,T}
    normal::SVector{2,T}
    id::Integer
    next::Integer
end
function InfiniteWall(sp::AbstractVector, ep::AbstractVector,
    n::AbstractVector, id::Integer, next::Integer)
    T = eltype(sp)
    n = normalize(n)
    d = dot(n, ep-sp)
    if abs(d) > 10eps(T)
        error("Normal vector is not actually normal to the wall")
    end
    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    return InfiniteWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), id, next)
end

"""
    FiniteWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing specular reflection during collision (immutable type).
Slower than [`InfiniteWall`](@ref), meant to be used for non-convex billiards.

Giving a `true` value to the field `isdoor` designates this obstacle to be a `Door`.
This is used in [`escapetime`](@ref) function. A `Door` is a obstacle of the
billiard that the particle can escape from, thus enabling calculations
of escape times.

### Fields:
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the
  particle *will come from before a collision* (pointing towards the inside of the
  billiard). The size of the vector is irrelevant
  since it is internally normalized.
* `isdoor::Bool` : Flag of whether this `FiniteWall` instance is a "Door".
  Defaults to `false`.
  * `id::Integer` : Obstacle id.
"""
struct FiniteWall{T<:AbstractFloat} <: Wall{T}
    sp::SVector{2,T}
    ep::SVector{2,T}
    normal::SVector{2,T}
    width::T
    center::SVector{2,T}
    isdoor::Bool
    id::Integer
    next::Integer
end
function FiniteWall(sp::AbstractVector, ep::AbstractVector,
    n::AbstractVector, isdoor::Bool = false, id::Integer, next::Integer)
    T = eltype(sp)
    n = normalize(n)
    d = dot(n, ep-sp)
    if abs(d) > 100eps(T)
        error("Normal vector is not actually normal to the wall: dot = $d")
    end
    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    w = norm(ep - sp)
    center = @. (ep+sp)/2
    return FiniteWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n),
    w, SVector{2,T}(center), isdoor, id, next)
end
FiniteWall(a, b, c, n::Integer, next::Integer) = FiniteWall(a, b, c, false, n, next)

isdoor(w) = w.isdoor

"""
    RandomWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing (uniformly) random reflection during collision (immutable type).
### Fields:
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the
  particle *is expected to come from* (pointing towards the inside of the
  billiard).
  * `id::Integer` : Obstacle id.
"""
struct RandomWall{T<:AbstractFloat} <: Wall{T}
    sp::SVector{2,T}
    ep::SVector{2,T}
    normal::SVector{2,T}
    id::Integer
    next::Integer
end
function RandomWall(sp::AbstractVector, ep::AbstractVector,
    n::AbstractVector, id::Integer, next::Integer)
    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    n = normalize(n)
    d = dot(n, ep-sp)
    if abs(d) > 10eps(T)
        error("Normal vector is not actually normal to the wall")
    end
    return RandomWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), id, next)
end

"""
    PeriodicWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle that imposes periodic boundary conditions upon collision (immutable type).
### Fields:
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the
  particle *will come from* (to the inside the billiard).
  The size of the vector is **important**!
  This vector is added to a particle's `pos` during collision. Therefore the
  size of the normal vector must be correctly associated with the size of the
  periodic cell.
  * `id::Integer` : Obstacle id.
"""
struct PeriodicWall{T<:AbstractFloat} <: Wall{T}
    sp::SVector{2,T}
    ep::SVector{2,T}
    normal::SVector{2,T}
    id::Integer
    next::Integer
end
function PeriodicWall(sp::AbstractVector, ep::AbstractVector,
    n::AbstractVector, id::Integer, next::Integer)
    T = eltype(sp)
    d = dot(n, ep-sp)
    if abs(d) > 10eps(T)
        error("Normal vector is not actually normal to the wall")
    end
    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    return PeriodicWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), id, next)
end


#pretty print:
show(io::IO, w::Wall{T}) where {T} = print(io, "Wall $(w.id) {$T}\n",
"start point: $(w.sp)\nend point: $(w.ep)\nnormal vector: $(w.normal)")

"""
    default_normal(sp, ep)
Return a vector to `v = ep - sp`, pointing to the left of `v`.
"""
function default_normal(sp, ep)
    T = eltype(sp)
    (x, y) = ep .- sp
    return SV{T}(-y, x)
end

for WT in (:InfiniteWall, :FiniteWall, :RandomWall)
    @eval $(WT)(sp, ep, id::Integer, next::Integer) = $(WT)(sp, ep, default_normal(sp, ep), id, next)
end


#######################################################################################
## Others: Ellipses
#######################################################################################

"""
    Ellipse{T<:AbstractFloat}  <: Obstacle{T}
Ellipse obstacle that also allows ray-splitting. The ellipse is always oriented
on the x and y axis (although you can make whichever you want the major one).
### Fields:
* `c::SVector{2,T}` : Center.
* `a::T` : x semi-axis.
* `b::T` : y semi-axis.
* `pflag::Bool` : Flag that keeps track of where the particle is currently
  propagating. `true` (default) is associated with being outside the ellipse.
  * `id::Integer` : Obstacle id.

The ellipse equation is given by
```math
\\left(\\frac{x - c[1]}{a} \\right)^2+ \\left(\\frac{y - c[2]}{b}\\right)^2 = 1
```
"""
mutable struct Ellipse{T<:AbstractFloat} <: Obstacle{T}
    c::SVector{2,T}
    a::T
    b::T
    pflag::Bool
    id::Integer
    next::Integer
    arc::T # arclength of one quadrant of the ellipse
end

"""
    ellipse_arclength(θ, e::Ellipse)
Return the arclength of the ellipse that
spans angle `θ` (in normal coordinates, not in the ellipse parameterization).
Expects `θ` to be in `[0, 2π]`.

After properly calculating the
```math
d=b\\,E\\bigl(\\tan^{-1}(a/b\\,\\tan(\\theta))\\,\\big|\\,1-(a/b)^2\\bigr)
```
"""
function ellipse_arclength(θ, e::Ellipse)
    n, θ = divrem(θ, π/2)
    a = e.a; b = e.b
    if a/b > 1.0
        θ = π/2 - θ
        return (n + 1.0)*e.arc - proper_ellipse_arclength(θ, b, a)
    else
        return n*e.arc + proper_ellipse_arclength(θ, a, b)
    end
end

proper_ellipse_arclength(θ, a, b)  = b*Elliptic.E(atan(a*tan(θ)/b), 1.0 - (a/b)^2)

export ellipse_arclength

function Ellipse(c::AbstractVector{T}, a, b, pflag = true,
                 id::Integer) where {T<:Real}
    S = T <: Integer ? Float64 : T
    return Ellipse{S}(SVector{2,S}(c),
            convert(S, abs(a)), convert(S, abs(b)), pflag, id, id,
            proper_ellipse_arclength(π/2, min(a, b), max(a, b))
        )
end
Ellipse{T}(args...) where {T} = Ellipse(args...)

#######################################################################################normalvec
## Normal vectors
#######################################################################################
"""
    normalvec(obst::Obstacle, position)
Return the vector normal to the obstacle's boundary at the given position (which is
assumed to be very close to the obstacle's boundary).
"""
@inline normalvec(wall::Wall, pos) = wall.normal
@inline normalvec(w::PeriodicWall, pos) = normalize(w.normal)
@inline normalvec(disk::Circular, pos) = normalize(pos - disk.c)
@inline normalvec(a::Antidot, pos) =
    a.pflag ? normalize(pos - a.c) : -normalize(pos - a.c)
@inline normalvec(d::Semicircle, pos) = normalize(d.c - pos)

@inline function normalvec(e::Ellipse{T}, pos) where {T}
    # from https://www.algebra.com/algebra/homework/
    # Quadratic-relations-and-conic-sections/Tangent-lines-to-an-ellipse.lesson
    x₀, y₀ = pos
    h, k = e.c
    s = e.pflag ? one(T) : -one(T)
    return s*normalize(SV((x₀-h)/(e.a*e.a), (y₀-k)/(e.b*e.b)))
end

#######################################################################################
## Distances
#######################################################################################
"""
    project_to_line(point, c, n)
Project given `point` to line that contains point `c` and has **normal vector** `n`.
"""
@inline function project_to_line(point, c, n)
    posdot = dot(c - point, n)
    intersection = point + posdot*n
end

"""
    distance(p::AbstractParticle, o::Obstacle)
Return the **signed** distance between particle `p` and obstacle `o`, based on
`p.pos`. Positive distance corresponds to the particle being on the *allowed* region
of the `Obstacle`. E.g. for a `Disk`, the distance is positive when the particle is
outside of the disk, negative otherwise.

    distance(p::AbstractParticle, bd::Billiard)
Return minimum `distance(p, obst)` for all `obst` in `bd`.
If the `distance(p, bd)` is negative and `bd` is convex,
this means that the particle is outside the billiard.

**WARNING** : `distance(p, bd)` may give negative values for non-convex
billiards, or billiards that are composed of several connected sub-billiards.

All `distance` functions can also be given a position (vector) instead of a particle.
"""
(distance(p::AbstractParticle{T}, obst::Obstacle{T})::T) where {T} =
distance(p.pos, obst)

@inline function distance(pos::AbstractVector{T}, w::Wall{T})::T where {T}
    v1 = pos - w.sp
    dot(v1, normalvec(w, pos))
end


@inline distance(pos::AbstractVector{T}, d::Circular{T}) where {T} =
    norm(pos - d.c) - d.r

@inline function distance(
    pos::AbstractVector{T}, a::Antidot{T})::T where {T}
    d = norm(pos - a.c) - a.r
    a.pflag ? d : -d
end

function distance(pos::AbstractVector{T}, s::Semicircle{T}) where {T}
    # Check on which half of circle is the particle
    v1 = pos - s.c
    nn = dot(v1, s.facedir)
    if nn ≤ 0 # I am "inside" semicircle
        return s.r - norm(pos - s.c)
    else # I am on the "other side"
        # end1 = SV(s.c[1] + s.r*s.facedir[2], s.c[2] - s.r*s.facedir[1])
        # end2 = SV(s.c[1] - s.r*s.facedir[2], s.c[2] + s.r*s.facedir[1])
        # return min(norm(pos - end1), norm(pos - end2))
        return one(T)
    end
end

function distance(pos::SV, e::Ellipse{T})::T where {T}
    d = ((pos[1] - e.c[1])/e.a)^2 + ((pos[2]-e.c[2])/e.b)^2 - 1.0
    e.pflag ? d : -d
end

# The entire functionality of `distance_init` is necessary only for
# FiniteWall !!!
distance_init(p::AbstractParticle, a::Obstacle) = distance_init(p.pos, a)
distance_init(pos::SVector, a::Obstacle) = distance(pos, a)

function distance_init(pos::SVector{2,T}, w::FiniteWall{T})::T where {T}

    n = normalvec(w, pos)
    posdot = dot(w.sp .- pos, n)
    if posdot ≥ 0 # I am behind wall
        intersection = project_to_line(pos, w.center, n)
        dfc = norm(intersection - w.center)
        if dfc > w.width/2
            return +1.0 # but not directly behind
        else
            return -1.0
        end
    end
    v1 = pos .- w.sp
    dot(v1, n)
end

####################################################
## Initial Conditions
####################################################
"""
    cellsize(bd)
Return the delimiters `xmin, ymin, xmax, ymax` of the given obstacle/billiard.

Used in `randominside()`, error checking and plotting.
"""
function cellsize(d::Circular{T}) where {T}
    xmin = ymin = T(Inf)
    xmax = ymax = T(-Inf)
    return xmin, ymin, xmax, ymax
end

function cellsize(w::Wall)
    xmin = min(w.sp[1], w.ep[1])
    xmax = max(w.sp[1], w.ep[1])
    ymin = min(w.sp[2], w.ep[2])
    ymax = max(w.sp[2], w.ep[2])
    return xmin, ymin, xmax, ymax
end

function cellsize(a::Antidot{T}) where {T}
    if a.pflag
        xmin = ymin = T(Inf)
        xmax = ymax = T(-Inf)
    else
        xmin, ymin = a.c .- a.r
        xmax, ymax = a.c .+ a.r
    end
    return xmin, ymin, xmax, ymax
end

function cellsize(e::Ellipse{T}) where {T}
    if e.pflag
        xmin = ymin = T(Inf)
        xmax = ymax = T(-Inf)
    else
        xmin = e.c[1] .- e.a; ymin = e.c[2] .- e.b
        xmax = e.c[1] .+ e.a; ymax = e.c[2] .+ e.b
    end
    return xmin, ymin, xmax, ymax
end

function cellsize(a::Semicircle{T}) where {T}
    xmin, ymin = a.c .- a.r
    xmax, ymax = a.c .+ a.r
    return xmin, ymin, xmax, ymax
end


####################################################
## Translations
####################################################
"""
    translate(obst::Obstacle, vector)
Create a copy of the given obstacle with its position
translated by `vector`.
"""
function translate end

for T in subtypes(Circular)
  @eval translate(d::$T, vec) = ($T)(d.c + vec, d.r, d.id, d.next)
end

for T in subtypes(Wall)
  @eval translate(w::$T, vec) = ($T)(w.sp + vec, w.ep + vec, w.normal, w.id, w.next)
end

translate(e::Ellipse, vec) = Ellipse(e.c + vec, e.a, e.b, e.id, e.next)


####################################################
## Previous obstacle
####################################################

function prev_obst(o::Obstacle{T}, bd::Billiard) where {T}
    for o2 in bd
        if o2.next == o.id
            return o2.next
        end
    end
    return 0
end