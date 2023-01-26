export reduce_boundary, add_collision

function add_collision(wall_hit::Obstacle{T}, new_wall::Obstacle{T}, vel::SV{T}) where {T}
    ending = new_wall.ep
    paramet_collision = dot(wall_hit.ep-wall_hit.sp, ending-wall_hit.sp)/(norm(wall_hit.ep-wall_hit.sp))^2
    moving_foward = dot(vel, wall_hit.ep - wall_hit.sp) > 0
    wall_hit_upside = dot(vel, wall_hit.normal) > 0
    incoming_upside = dot(vel, new_wall.normal) < 0

    return (wall_hit.id, new_wall.id, paramet_collision, moving_foward, wall_hit_upside, incoming_upside)
end

function search_next(searchparams::NamedTuple, collisions::DataFrame, bd::Billiard{T}) where {T}
    if searchparams[:moving_foward]
        selection = collisions[(collisions.wall_hit_id .== searchparams[:wall_hit]) .& (collisions.wall_hit_upside .== searchparams[:wall_hit_upside]).& (collisions.collision_point .> searchparams[:collision_point]), :]
    else
        selection = collisions[(collisions.wall_hit_id .== searchparams[:wall_hit]) .& (collisions.wall_hit_upside .== searchparams[:wall_hit_upside]).& (collisions.collision_point .< searchparams[:collision_point]), :]
    end

    for obst in bd
        if obst.id == searchparams[:wall_hit]
            current_wall = obst
            break
        end
    end

    if isempty(selection)
        if searchparams[:moving_foward]
            new_col = collisions[(collisions.incoming_id .== searchparams[:wall_hit]), :]
            # Special case: outer wall hit
            if new_col.wall_hit_id[1] == 0
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
            new_col = collisions[collisions.incoming_id .== prev_obst(current_wall, bd), :]
            # Special case: outer wall hit
            if new_col.wall_hit_id[1] == 0
                next_wall = prev_obst(current_wall)
                new_collision_point = 1.0
                travel_dir = false
                next_side = true
                return (wall_hit = next_wall, collision_point = new_collision_point, moving_foward = travel_dir, wall_hit_upside = next_side)
            end

            # Special case: first wall moving backwards
            if current_wall.id == current_wall.next
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

        for obst in bd
            if obst.id == new_col.incoming_id[1]
                incoming_obst = obst
            end
        end

        if new_col.moving_foward[1] == searchparams[:moving_foward]
            next_wall = incoming_obst.id
            new_collision_point = 1.0
            travel_dir = false
            next_side = new_col.incoming_upside[1]
        else
            next_wall = incoming_obst.next
            new_collision_point = 0.0
            travel_dir = true
            next_side = new_col.incoming_upside[1]
        end
    end

    return (wall_hit = next_wall, collision_point = new_collision_point, moving_foward = travel_dir, wall_hit_upside = next_side)
end

function reduce_boundary!(bd_active::Billiard, collisions::DataFrame, current_searchparams::NamedTuple)
    new_shape_ids = Set()
    bd_new = Obstacle[]
    current_wall = last(bd_active).id
    push!(new_shape_ids, current_wall)
    
    counter = 0
    while true
        push!(new_shape_ids, current_searchparams[:wall_hit])
        current_searchparams = search_next(current_searchparams, collisions, sides)
        if current_searchparams[:wall_hit] == current_wall || counter > 20
            break
        end
        counter += 1
    end
    # collisions_new = filter(row -> (row.wall_hit_id ∈ new_shape_ids || row.incoming_id ∈ new_shape_ids || row.incoming_id ∈ (new_shape_ids .+ 1)), collisions)

    for obst in bd_active
        if parse(Int64, obst.name) in new_shape_ids
            push!(bd_new, obst)
        end
    end
    return Billiard(bd_new)
end