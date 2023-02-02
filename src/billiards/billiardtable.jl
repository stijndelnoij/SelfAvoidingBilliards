export Billiard, randominside, randominside_xyφ, prev_obst
#######################################################################################
## Billiard Table
#######################################################################################
mutable struct Billiard{T, O<:Array, Y}
    obstacles::O
    peridx::Y
end

@inline isperiodic(::Billiard{T,<:Array,Vector{Int}}) where {T} = true
@inline isperiodic(::Billiard{T,<:Array,Nothing}) where {T} = false
@inline isperiodic(i::Int, ::Billiard{T,O,Nothing}) where {T,O} = false
@inline isperiodic(i::Int, bd::Billiard{T,<:Array,Vector{Int}}) where {T} = i ∈ bd.peridx

#pretty print:
_get_name(o::Obstacle) = :name ∈ fieldnames(typeof(o)) ? o.id : string(typeof(o))
function Base.show(io::IO, bd::Billiard{T,BT}) where {T, BT}
    D = length(bd)
    s = "Billiard{$T} with $D obstacles:"
    bd.peridx != nothing && (s = "Periodic "*s)
    for i in 1:min(10, D)
        s*="\n  $(_get_name(bd[i]))"
    end
    print(io, s)
    length(bd) > 10 && print(io, "\n  ...")
end



"""
    Billiard(obstacles...)
Construct a `Billiard` from given `obstacles` (tuple, vector, varargs).

For functions like [`boundarymap`](@ref),
it is expected (if possible) that the obstacles of the billiard are sorted,
such that the arc-coordinate `ξ` around the billiard is increasing counter-clockwise.

`ξ` is measured as:
* the distance from start point to end point in `Wall`s
* the arc length measured counterclockwise from the open face in `Semicircle`s
* the arc length measured counterclockwise from the rightmost point
  in `Circular`s
"""
function Billiard(bd::Union{AbstractVector, Tuple})

    T = eltype(bd[1])
    D = length(bd)
    # Assert that all elements of `bd` are of same type:
    for i in 2:D
        eltype(bd[i]) != T && throw(ArgumentError(
        "All obstacles of the billiard must have same type of
        numbers. Found $T and $(eltype(bd[i])) instead."
        ))
    end
    tup = Array(bd)
    peridx = findall(x -> typeof(x) <: PeriodicWall, tup)
    if isodd(length(peridx))
        throw(ArgumentError(
        "A billiard can only have an even number of `PeriodicWall`s, "*
        "since they have to come in pairs."
        ))
    end

    if !any(x -> x isa PeriodicWall, bd)
        return Billiard{T, typeof(tup), Nothing}(tup, nothing)
    else
        return Billiard{T, typeof(tup), Vector{Int}}(tup, peridx)
    end
end

function Billiard(bd::Vararg{Obstacle})
    tup = [bd...,]
    return Billiard(tup)
end


Base.getindex(bd::Billiard, i) = bd.obstacles[i]
Base.lastindex(bd::Billiard) = length(bd.obstacles)
Base.firstindex(bd::Billiard) = 1
# Iteration:
Base.iterate(bd::Billiard) = iterate(bd.obstacles)
Base.iterate(bd::Billiard, state) = iterate(bd.obstacles, state)
Base.length(bd::Billiard) = length(bd.obstacles) 
# Base.length(::Billiard{T, L}) where {T, L} = L
Base.eltype(::Billiard{T}) where {T} = T

#######################################################################################
## Distances
#######################################################################################
for f in (:distance, :distance_init)
    @eval $(f)(p::AbstractParticle, bd::Billiard) = $(f)(p.pos, bd.obstacles)
    @eval $(f)(pos::SV{T}, bd::Billiard) where {T} = $(f)(pos, bd.obstacles)
end

for f in (:distance, :distance_init)
    @eval begin
        function ($f)(p::SV{T}, bd::Array)::T where {T}
            dmin::T = T(Inf)
            for obst in bd
                d::T = distance(p, obst)
                d < dmin && (dmin = d)
            end
            return dmin
        end
    end
end

#######################################################################################
## randominside
#######################################################################################
function cellsize(
    bd::Union{Vector{<:Obstacle{T}}, Billiard{T}}) where {T<:AbstractFloat}

    xmin::T = ymin::T = T(Inf)
    xmax::T = ymax::T = T(-Inf)
    for obst ∈ bd
        xs::T, ys::T, xm::T, ym::T = cellsize(obst)
        xmin = xmin > xs ? xs : xmin
        ymin = ymin > ys ? ys : ymin
        xmax = xmax < xm ? xm : xmax
        ymax = ymax < ym ? ym : ymax
    end
    return xmin, ymin, xmax, ymax
end

"""
    randominside(bd::Billiard) → p
Return a particle `p` with random allowed initial conditions inside the given
billiard.

**WARNING** : `randominside` works for any **convex** billiard but it does
not work for non-convex billiards. (this is because it uses `distance`)
"""
randominside(bd::Billiard) = Particle(_randominside(bd)...)

"""
    randominside_xyφ(bd::Billiard) → x, y, φ
Same as [`randominside`](@ref) but only returns position and direction.
"""
function randominside_xyφ(bd::Billiard{T}) where {T<:AbstractFloat}
    xmin::T, ymin::T, xmax::T, ymax::T = cellsize(bd)
    xp = T(rand())*(xmax-xmin) + xmin
    yp = T(rand())*(ymax-ymin) + ymin
    pos = SV{T}(xp, yp)
    dist = distance_init(pos, bd)
    while dist <= sqrt(eps(T))
        xp = T(rand())*(xmax-xmin) + xmin
        yp = T(rand())*(ymax-ymin) + ymin
        pos = SV{T}(xp, yp)
        dist = distance_init(pos, bd)
    end
    return pos[1], pos[2], T(2π*rand())
end

function _randominside(bd::Billiard{T}) where {T<:AbstractFloat}
    x, y, φ = randominside_xyφ(bd)
    pos = SV{T}(x, y)
    vel = SV{T}(sincos(φ)...)
    return pos, vel, zero(SV{T}) # this zero is the `current_cell`
end

####################################################
## Previous obstacle
####################################################

function prev_obst(o::Obstacle{T}, bd::Billiard{T}) where {T}
    prev = filter(e -> e.next == o.id, bd.obstacles)
    if isempty(prev)
        return 0
    else
        return first(prev).id
    end
end