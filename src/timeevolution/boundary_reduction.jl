export reduce_boundary!, add_collision, CollisionList

"""
    add_collision(wall_hit::Obstacle, new_wall::Obstacle, vel) -> new_collision
Returns a tuple containing a row for the `collisions` dataframe.

The tuple contains:
* index of the wall which the particle hit
* index of the newly placed wall
* parametric position on the wall where the particle has collided (between 0 and 1)
* boolean wether the particle will move foward with respect to the wall after the collision
* boolean wether the particle hit the "upside" of the wall
* boolean wether the upside of the new wall is facing the collision
"""
function add_collision(wall_hit::Obstacle{T}, new_wall::Obstacle{T}, vel::SV{T}) where {T}
    ending = new_wall.ep
    if wall_hit isa Wall
        paramet_collision = dot(wall_hit.ep-wall_hit.sp, ending-wall_hit.sp)/(norm(wall_hit.ep-wall_hit.sp))^2
        moving_foward = dot(vel, wall_hit.ep - wall_hit.sp) > 0
    elseif wall_hit isa Ellipse || wall_hit isa Circular
        paramet_collision = (atan(ending[2] - wall_hit.c[2], ending[1]- wall_hit.c[1])+π) / (2*π)
        moving_foward = dot(vel, [-ending[2], ending[1]]) > 0
    # If not wall/circle, no support for boundary reduction
    else
        paramet_collision = 0.
        moving_foward = true
    end
    wall_hit_upside = dot(vel, normalvec(wall_hit, ending)) > 0
    incoming_upside = dot(vel, normalvec(new_wall, ending)) < 0

    return (wall_hit.id, new_wall.id, paramet_collision, moving_foward, wall_hit_upside, incoming_upside)
end

# Searches next wall connected
function search_next(searchparams::NamedTuple, collisions::DataFrame, bd::Billiard{T}) where {T}
    if searchparams[:moving_foward]
        selection = collisions[(collisions.wall_hit_id .== searchparams[:wall_hit]) .& (collisions.wall_hit_upside .== searchparams[:wall_hit_upside]) .& (collisions.collision_point .> searchparams[:collision_point]), :]
    else
        selection = collisions[(collisions.wall_hit_id .== searchparams[:wall_hit]) .& (collisions.wall_hit_upside .== searchparams[:wall_hit_upside]).& (collisions.collision_point .< searchparams[:collision_point]), :]
    end

    # Return error code if wall searching is not in the active billiard table
    if isempty(filter(e -> e.id == searchparams[:wall_hit], bd.obstacles))
        return (wall_hit = 0, collision_point = 0., moving_foward = true, wall_hit_upside = true)
    end

    current_wall = first(filter(e -> e.id == searchparams[:wall_hit], bd.obstacles))

    if isempty(selection)
        if searchparams[:moving_foward]
            new_col = collisions[(collisions.incoming_id .== searchparams[:wall_hit]), :]
            # Special case: outer wall hit
            if !(current_wall isa Tail)
                next_wall = current_wall.next
                new_collision_point = 0.0
                travel_dir = true
                next_side = true
                return (wall_hit = next_wall, collision_point = new_collision_point, moving_foward = travel_dir, wall_hit_upside = next_side)
            end

            if new_col.incoming_upside[1] == searchparams[:wall_hit_upside]
                next_wall = new_col.wall_hit_id[1]
                new_collision_point = new_col.collision_point[1]
                travel_dir = !new_col.moving_foward[1]
                next_side = new_col.wall_hit_upside[1]
            else
                next_wall = current_wall.next
                new_collision_point = 0.0
                travel_dir = true
                next_side = !new_col.incoming_upside[1]
            end
        else
            new_col = collisions[collisions.incoming_id .== current_wall.id - 1, :]
            # Special case: outer wall hit
            if !(current_wall isa Tail)
                next_wall = prev_obst(current_wall, bd)
                new_collision_point = 1.0
                travel_dir = false
                next_side = true
                return (wall_hit = next_wall, collision_point = new_collision_point, moving_foward = travel_dir, wall_hit_upside = next_side)
            end

            # Special case: first wall moving backwards
            if isempty(new_col)
                next_wall = current_wall.id
                new_collision_point = 0.0
                travel_dir = true
                next_side = !searchparams[:wall_hit_upside]
                return (wall_hit = next_wall, collision_point = new_collision_point, moving_foward = travel_dir, wall_hit_upside = next_side)
            end

            

            if new_col.incoming_upside[1] == searchparams[:wall_hit_upside]
                next_wall = new_col.wall_hit_id[1]
                new_collision_point = new_col.collision_point[1]
                travel_dir = new_col.moving_foward[1]
                next_side = new_col.wall_hit_upside[1]
            else
                next_wall = new_col.incoming_id[1]
                new_collision_point = 1.0
                travel_dir = false
                next_side = searchparams[:wall_hit_upside]
            end
        end
    else
        if searchparams[:moving_foward]
            new_col = selection[findmin(selection.collision_point)[2], :]
        else
            new_col = selection[findmax(selection.collision_point)[2], :]
        end

        # This below only works for Tail

        if new_col.moving_foward[1] == searchparams[:moving_foward]
            next_wall = new_col.incoming_id[1]
            new_collision_point = 1.0
            travel_dir = false
            next_side = new_col.incoming_upside[1]
        else
            next_wall = new_col.incoming_id[1] + 1
            new_collision_point = 0.0
            travel_dir = true
            next_side = new_col.incoming_upside[1]
        end
    end

    return (wall_hit = next_wall, collision_point = new_collision_point, moving_foward = travel_dir, wall_hit_upside = next_side)
end

"""
    reduce_boundary!(bd::Billiard, collisions::DataFrame) -> area
Reduces the current table to the smallest nessecary size by checking which walls are connected to the current particle's position.
Also returns the vertices of the enclosing polygon.

Note that this function will only work if all the obstacles are (segments of) lines. Please do not use when circular/elliptical obstacles are involved!
"""
function reduce_boundary!(bd_active::Billiard, collisions::DataFrame)
    area = 0.0
    vertices = Vector{Vector{Float64}}()
    last_collision = last(collisions)
    current_searchparams = (wall_hit = last_collision[1], collision_point = last_collision[3], moving_foward = last_collision[4], wall_hit_upside = last_collision[5])
    new_shape_ids = Set()
    current_wall = last(bd_active).id
    push!(new_shape_ids, current_wall)
    obst = filter(e -> e.id == current_searchparams[:wall_hit], bd_active.obstacles)[1]
    vertex = (obst.ep - obst.sp) * current_searchparams[:collision_point] + obst.sp
    push!(vertices, vertex)
    
    counter = 0

    prev_searchparams = current_searchparams
    while true
        push!(new_shape_ids, current_searchparams[:wall_hit])
        current_searchparams = search_next(current_searchparams, collisions, bd_active)
        
        # if current wall not in the active table, then there was a problem like perpendicular collision or invalid backpropagation
        if current_searchparams[:wall_hit] ∈ [o.id for o in bd_active.obstacles]
            obst = filter(e -> e.id == current_searchparams[:wall_hit], bd_active.obstacles)[1]
        else
            return -1.
        end
        vertex = (obst.ep - obst.sp) * current_searchparams[:collision_point] + obst.sp
        push!(vertices, vertex)
        if current_searchparams[:wall_hit] == current_wall || counter > 2*length(bd_active) || current_searchparams[:wall_hit] == 0
            break
        end
        counter += 1
        prev_searchparams = current_searchparams
    end

    if current_searchparams[:wall_hit] != 0
        filter!(e -> e.id ∈ new_shape_ids, bd_active.obstacles)
    end
    return vertices
end

"""
    CollisionList()
    Returns empty collisions dataframe, use this when initializing a SAB run.
"""
function CollisionList()
    return DataFrame(wall_hit_id = Integer[], incoming_id = Integer[], collision_point = AbstractFloat[], moving_foward = Bool[], wall_hit_upside = Bool[], incoming_upside = Bool[])
end