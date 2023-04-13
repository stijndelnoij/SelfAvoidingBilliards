using StaticArrays

export billiard_rectangle, billiard_sinai, billiard_polygon, billiard_lorentz, 
billiard_hexagonal_sinai, billiard_bunimovich,
billiard_stadium, billiard_mushroom, billiard_logo, billiard_vertices,
polygon_vertices, billiard_deformed_triangle

####################################################
## Famous/Standard Billiards
####################################################
"""
    billiard_rectangle(x=1.0, y=1.0; setting = "standard")
Return a vector of obstacles that defines a rectangle billiard of size (`x`, `y`).

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type,
  enforcing periodicity at the boundaries
* "random" : The velocity is randomized upon collision.
* "ray-splitting" : All obstacles in the billiard allow for ray-splitting.
"""
function billiard_rectangle(x′ = 1.0, y′ = 1.0;
    x = x′, y = y′, setting::String = "standard")

    x = convert(AbstractFloat, x)
    x, y = promote(x,y)
    o = typeof(x)(0.0)
    if setting == "standard"
        sp = [o,y]; ep = [o, o]; n = [x,o]
        leftw = InfiniteWall(sp, ep, n, 1, 2)
        sp = [x,o]; ep = [x, y]; n = [-x,o]
        rightw = InfiniteWall(sp, ep, n, 3, 4)
        sp = [x,y]; ep = [o, y]; n = [o,-y]
        topw = InfiniteWall(sp, ep, n, 2, 3)
        sp = [o,o]; ep = [x, o]; n = [o,y]
        botw = InfiniteWall(sp, ep, n, 4, 1)
    elseif setting == "periodic"
        sp = [o,y]; ep = [o, o]; n = [x,o]
        leftw = PeriodicWall(sp, ep, n, 1, 2)
        sp = [x,o]; ep = [x, y]; n = [-x,o]
        rightw = PeriodicWall(sp, ep, n, n, 3, 4)
        sp = [x,y]; ep = [o, y]; n = [o,-y]
        topw = PeriodicWall(sp, ep, n, 2, 3)
        sp = [o,o]; ep = [x, o]; n = [o,y]
        botw = PeriodicWall(sp, ep, n, 4, 1)
    elseif setting == "random"
        sp = [o,y]; ep = [o, o]; n = [x,o]
        leftw = RandomWall(sp, ep, n, 1, 2)
        sp = [x,o]; ep = [x, y]; n = [-x,o]
        rightw = RandomWall(sp, ep, n, n, 3, 4)
        sp = [x,y]; ep = [o, y]; n = [o,-y]
        topw = RandomWall(sp, ep, n, 2, 3)
        sp = [o,o]; ep = [x, o]; n = [o,y]
        botw = RandomWall(sp, ep, n, 4, 1)
    else
        throw(ArgumentError("The given setting=$setting is unknown."))
    end
    return Billiard(botw, rightw, topw, leftw)
end

"""
    billiard_sinai(r=0.25, x=1.0, y=1.0; setting = "standard")
Return a vector of obstacles that defines a Sinai billiard of size (`x`, `y`) with
a disk in its center, of radius `r`.

In the periodic case, the system is also known as "Lorentz Gas".

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type,
  enforcing periodicity at the boundaries
* "random" : The velocity is randomized upon collision.
"""
function billiard_sinai(r′ = 0.25, x′ = 1.0, y′ = 1.0;
    r = r′, x = x′, y = y′, setting = "standard")

    if (setting == "periodic") && (r>=x/2 || r>=y/2)
        es = "Disk radius too big for a periodic Sinai billiard.\n"
        es*= "Obstacles must not overlap with `PeriodicWall`s."
        error(es)
    end
    r, x, y = promote(r,x,y)
    bdr = billiard_rectangle(x, y; setting = setting)

    c = [x/2, y/2]
    if setting == "random"
        centerdisk = RandomDisk(c, r, 5)
    else
        centerdisk = Disk(c, r, 5)
    end

    return Billiard(centerdisk, bdr...)
end

"""
    billiard_lorentz(r=0.25, x=1.0, y=1.0)
Alias for `billiard_sinai(r,x,y; setting = "periodic")`.
"""
billiard_lorentz(r=0.25, x=1.0, y=1.0) = billiard_sinai(r,x,y; setting = "periodic")

"""
    billiard_polygon(n::Int, R, center = [0,0]; setting = "standard")
Return a vector of obstacles that defines a regular-polygonal billiard
with `n` sides, radius `r` and given `center`.

Note: `R` denotes the so-called outer radius, not the inner one.

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type, enforcing periodicity
  at the boundaries. Only available for `n=4` or `n=6`.
* "random" : The velocity is randomized upon collision.
"""
function billiard_polygon(sides′::Int, r′::Real, center′ = [0,0];
    sides::Int = sides′, r = r′, center = center′, setting = "standard")

    S = typeof(convert(AbstractFloat, r))
    bd = Obstacle{S}[]
    verteces = polygon_vertices(r, sides, center)

    if setting == "standard"
        T = InfiniteWall
    elseif setting == "periodic"
        if sides != 4 && sides != 6
            error("Polygonal and periodic billiard can exist only for `n=4` or `n=6`")
        end
        T = PeriodicWall
        inr = S(r*cos(π/sides))
    elseif setting == "random"
        T = RandomWall
    end

    for i in eachindex(verteces)
        starting = verteces[i]
        ending = verteces[mod1(i+1, sides)]
        # Normal vector must look at where the particle is coming from
        w = ending - starting
        if setting == "periodic"
            normal = 2inr*normalize([-w[2], w[1]])
            wall = T(starting, ending, normal, i, mod1(i+1, sides))
        else
            normal = [-w[2], w[1]]
            wall = T(starting, ending, normal, i, mod1(i+1, sides))
        end
        push!(bd, wall)
    end
    return Billiard(bd)
end

"""
    Return deformed equallateral triangle billiard. ϵ is the deformation vector. 
"""

function billiard_deformed_triangle(ϵ, r::Real, center = [0,0])
    S = typeof(convert(AbstractFloat, r))
    bd = Obstacle{S}[]
    sides = 3
    verteces = polygon_vertices(r, sides, center, -π/6)
    verteces[1] = verteces[1] + ϵ

    T = InfiniteWall

    for i in eachindex(verteces)
        starting = verteces[i]
        ending = verteces[mod1(i+1, sides)]
        # Normal vector must look at where the particle is coming from
        w = ending - starting
        normal = [-w[2], w[1]]
        wall = T(starting, ending, normal, i, mod1(i+1, sides))
        push!(bd, wall)
    end
    return Billiard(bd)
end

"""
    polygon_vertices(r, sides, center = [0, 0], θ=0) -> v
Return the vertices that define a regular polygon with `sides` sides and
radius `r`, centered at `center`.

Optionally rotate by `θ` degrees.
"""
function polygon_vertices(r, sides, center = [0, 0.0], θ = 0)
    S = typeof(convert(AbstractFloat, r))
    [S[r*cos(θ + 2π*i/sides), r*sin(θ + 2π*i/sides)] .+ center for i in 1:sides]
end




"""
    billiard_hexagonal_sinai(r, R, center = [0,0]; setting = "standard")
Create a sinai-like billiard, which is a hexagon of outer radius `R`, containing
at its center (given by `center`) a disk of radius `r`. The `setting` keyword
is passed to `billiard_polygon`.
"""
function billiard_hexagonal_sinai(r′::Real = 0.5, R′::Real = 1.0, center′ = [0,0];
    r = r′, R = R′, center = center′, setting = "standard")
    r, R = promote(r, R)
    T = typeof(r); center = T[center...]
    bdr = billiard_polygon(6, R, center; setting = setting)
    DT = setting == "random" ? RandomDisk : Disk
    return Billiard(Disk(center, r, 7), bdr...)
end





"""
    billiard_bunimovich(l=1.0, w=1.0)
Return a vector of `Obstacle`s that define a Buminovich billiard, also called a
stadium. The length is considered *without* the attached semicircles, meaning that the
full length of the billiard is `l + w`. The left and right edges of the stadium
are [`Semicircle`](@ref)s.

`billiard_stadium` is an alias of `billiard_bunimovich`.
"""
function billiard_bunimovich(l′=1.0, w′=1.0; l = l′, w = w′)

    l = convert(AbstractFloat, l)
    l, w = promote(l,w)
    o = typeof(l)(0.0)
    bw = InfiniteWall([o,o], [l,o], [o,  w], 4, 1)
    tw = InfiniteWall([l,w], [o,w], [o, -w], 2, 3)
    leftc = Semicircle([o, w/2], w/2, [l, o], 1, 2)
    rightc = Semicircle([l, w/2], w/2, [-l, o], 3, 4)
    return Billiard(bw, rightc, tw, leftc)
end

billiard_stadium = billiard_bunimovich


"""
    billiard_vertices(v, type = FiniteWall)
Construct a polygon billiard that connects the given vertices `v`
(vector of 2-vectors). The vertices should construct a billiard in a
counter-clockwise orientation (i.e. the normal vector always points to the
left of `v[i+1] - v[i]`.).

`type` decides what kind of walls to use.
Ths function assumes that the first entry of `v` should be connected with the
last.
"""
function billiard_vertices(vertices, type = FiniteWall)
    @assert vertices[end] != vertices[1]
    obs = type[]
    for i in 1:length(vertices)-1
        sp, ep = vertices[i], vertices[i+1]
        push!(obs, type(sp, ep, i, i+1))
    end
    sp, ep = vertices[end], vertices[1]
    push!(obs, type(sp, ep, length(vertices), 1))
    return Billiard(obs)
end
