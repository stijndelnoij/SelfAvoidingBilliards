export reduce_boundary!, add_collision, CollisionList

function add_collision(wall_hit::Obstacle{T}, new_wall::Obstacle{T}, vel::SV{T}) where {T}
    ending = new_wall.ep
    paramet_collision = dot(wall_hit.ep-wall_hit.sp, ending-wall_hit.sp)/(norm(wall_hit.ep-wall_hit.sp))^2
    moving_foward = dot(vel, wall_hit.ep - wall_hit.sp) > 0
    wall_hit_upside = dot(vel, wall_hit.normal) > 0
    incoming_upside = dot(vel, new_wall.normal) < 0

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

function reduce_boundary!(bd_active::Billiard, collisions::DataFrame)
    last_collision = last(collisions)
    current_searchparams = (wall_hit = last_collision[1], collision_point = last_collision[3], moving_foward = last_collision[4], wall_hit_upside = last_collision[5])
    new_shape_ids = Set()
    current_wall = last(bd_active).id
    push!(new_shape_ids, current_wall)
    
    counter = 0
    while true
        push!(new_shape_ids, current_searchparams[:wall_hit])
        current_searchparams = search_next(current_searchparams, collisions, bd_active)
        if current_searchparams[:wall_hit] == current_wall || counter > 2*length(bd_active) || current_searchparams[:wall_hit] == 0
            break
        end
        counter += 1
    end
    if current_searchparams[:wall_hit] != 0
        filter!(e -> e.id ∈ new_shape_ids, bd_active.obstacles)
    end
end

"""
    CollisionList gives empty collisionlist dataframe
"""
function CollisionList()
    return DataFrame(wall_hit_id = Integer[], incoming_id = Integer[], collision_point = AbstractFloat[], moving_foward = Bool[], wall_hit_upside = Bool[], incoming_upside = Bool[])
end