using LinearAlgebra
using StaticArrays
using Elliptic
using Statistics
using DataFrames
using CSV
using Printf
include("DynamicalBilliards.jl-master/src/DynamicalBilliards.jl")
using .DynamicalBilliards2

const accuracy = 1e-4
const timestep = 0.001*sqrt(3)
const max_length = 6*sqrt(3)
const N_length_bins = 60
const d_length = max_length/N_length_bins

function search_next(searchparams::NamedTuple, collisions::DataFrame, sides::Int64)
    if searchparams[:moving_foward]
        selection = collisions[(collisions.wall_hit_id .== searchparams[:wall_hit]) .& (collisions.wall_hit_upside .== searchparams[:wall_hit_upside]).& (collisions.collision_point .> searchparams[:collision_point]), :]
    else
        selection = collisions[(collisions.wall_hit_id .== searchparams[:wall_hit]) .& (collisions.wall_hit_upside .== searchparams[:wall_hit_upside]).& (collisions.collision_point .< searchparams[:collision_point]), :]
    end

    if isempty(selection)
        if searchparams[:moving_foward]
            # Special case: outer wall hit
            if searchparams[:wall_hit] <= sides
                next_wall = mod1(searchparams[:wall_hit] + 1 , 3)
                new_collision_point = 0.0
                travel_dir = true
                next_side = true
                return (wall_hit = next_wall, collision_point = new_collision_point, moving_foward = travel_dir, wall_hit_upside = next_side)
            end

            new_col = collisions[(collisions.incoming_id .== searchparams[:wall_hit]), :]
            
            if new_col.incoming_upside[1] == searchparams[:wall_hit_upside]
                next_wall = new_col.wall_hit_id[1]
                new_collision_point = new_col.collision_point[1]
                travel_dir = !new_col.moving_foward[1]
                next_side = new_col.wall_hit_upside[1]
            else
                next_wall = new_col.incoming_id[1] + 1 # Could get larger than current wall !?
                new_collision_point = 0.0
                travel_dir = true
                next_side = !new_col.incoming_upside[1]
            end
        else
            # Special case: outer wall hit
            if searchparams[:wall_hit] <= sides
                next_wall = mod1(searchparams[:wall_hit] - 1 , 3)
                new_collision_point = 1.0
                travel_dir = false
                next_side = true
                return (wall_hit = next_wall, collision_point = new_collision_point, moving_foward = travel_dir, wall_hit_upside = next_side)
            end

            # Special case: first wall moving backwards
            if searchparams[:wall_hit] == sides + 1
                next_wall = sides + 1
                new_collision_point = 0.0
                travel_dir = true
                next_side = !searchparams[:wall_hit_upside]
                return (wall_hit = next_wall, collision_point = new_collision_point, moving_foward = travel_dir, wall_hit_upside = next_side)
            end

            new_col = collisions[collisions.incoming_id .== searchparams[:wall_hit] - 1, :]

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

        if new_col.moving_foward[1] == searchparams[:moving_foward]
            next_wall = new_col.incoming_id[1]
            new_collision_point = 1.0
            travel_dir = false
            next_side = new_col.incoming_upside[1]
        else
            next_wall = new_col.incoming_id + 1
            new_collision_point = 0.0
            travel_dir = true
            next_side = new_col.incoming_upside[1]
        end
    end

    return (wall_hit = next_wall, collision_point = new_collision_point, moving_foward = travel_dir, wall_hit_upside = next_side)
end

function reduce_boundary(bd_active::Billiard, iteration::Int64, collisions::DataFrame, current_searchparams::NamedTuple)
    new_shape_ids = Set()
    bd_new = Obstacle[]
    open_state = false
    push!(new_shape_ids, iteration)
    
    counter = 0
    while true
        push!(new_shape_ids, current_searchparams[:wall_hit])
        current_searchparams = search_next(current_searchparams, collisions, sides)
        if current_searchparams[:wall_hit] == sides + 1 && sides + 1 ∈ new_shape_ids
            open_state = true
        end
        if current_searchparams[:wall_hit] == iteration || counter > 20
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
    return Billiard(bd_new), open_state#, collisions_new
end


function add_collision(wall_hit::Obstacle, new_wall::Obstacle, vel::SArray{Tuple{2},Float64,1,2})
    ending = new_wall.ep
    paramet_collision = dot(wall_hit.ep-wall_hit.sp, ending-wall_hit.sp)/(norm(wall_hit.ep-wall_hit.sp))^2
    moving_foward = dot(vel, wall_hit.ep - wall_hit.sp) > 0
    wall_hit_upside = dot(vel, wall_hit.normal) > 0
    incoming_upside = dot(vel, new_wall.normal) < 0

    return (parse(Int64, wall_hit.name), parse(Int64, new_wall.name), paramet_collision, moving_foward, wall_hit_upside, incoming_upside)
end


function dist_from_start(wall_hit::Obstacle, coord::AbstractVector)
    if isa(wall_hit, SelfAvoidingBilliard)
        t0 = wall_hit.length
        dt = norm(coord - wall_hit.sp)
        return t0+dt
    else
        return -1.
    end
end


function stop_run(t::Float64, id_hit::Int64, α::Float64, bd_size::Int64)
    if t < accuracy
        return true, 0
    elseif id_hit == 0
        return true, 1
    elseif α < 1e-6
        return true, 2
    elseif bd_size > 1000
        return true, 3
    else
        return false, 0
    end
end


function end_pos(bd_active::Billiard, p::Particle, savepath::Bool = false)
    k = 1
    w = Array{Float64,1}(undef,2)
    # dt = Float64[]
    # angles = Float64[]
    # pathlengths = Float64[]

    error = Int64
    α = Float64
    total_t = 0.0
    new_searchparams = NamedTuple()
    collisions = DataFrame(wall_hit_id = Int64[], incoming_id = Int64[], collision_point = Float64[], moving_foward = Bool[], wall_hit_upside = Bool[], incoming_upside = Bool[])

    local start = p.pos
    local ending = Array{Float64,1}(undef,2)
    if savepath
        open("data/path.txt","w") do io
            println(io,start[1], " ", start[2])
        end
    end

    while true
        start = p.pos
        i, t, ending, vel=bounce!(p,bd_active)            
        
        if i != 0
            wall_hit = bd_active[i]
            # t_hit = dist_from_start(wall_hit, ending)
            # if t_hit >= 0
            #     push!(dt, t+total_t-t_hit)
            # end

            dotpr = abs(dot(vel, wall_hit.normal))/(norm(vel)*norm(wall_hit.normal))
            if dotpr <= 1.
                α = acos(dotpr)
            else
                α = acos(1.)
            end
        end

        w = ending - start
        normal = [-w[2], w[1]]

        # error codes: 0 no error, 1 no collision detected, 2 perpendicular collision, 3 too large bdy
        final_pos_reached, error = stop_run(t, i, α, k)

        if final_pos_reached
            break
        end

        # Save positions based on timeinterval
        # if total_t % timestep + t > timestep
        #     new_pos_add = [start + a*w/norm(w) for a in (timestep - total_t % timestep):timestep:t]
        #     append!(positions_t, new_pos_add)
        # end
        
        # push!(angles, α)
        # push!(pathlengths, t)

        if savepath
            open("data/path.txt","a") do io
                println(io,ending[1], " ", ending[2])
            end
        end

        # Adding new wall
        new_id_wall = string(k+sides)
        try
            new_wall = SelfAvoidingBilliard(start, ending, normal, total_t, new_id_wall)
        catch e
            final_pos_reached = true
            error = 2
            break
        end
        new_wall = SelfAvoidingBilliard(start, ending, normal, total_t, new_id_wall)
        # push!(bd_active, new_wall)
        bd_active = Billiard((bd_active..., new_wall))

        new_collision = add_collision(wall_hit, new_wall, vel)
        push!(collisions, new_collision)
        total_t += t

        if k % 20== 0
            new_searchparams = (wall_hit = new_collision[1], collision_point = new_collision[3], moving_foward = new_collision[4], wall_hit_upside = new_collision[5])
            bd_active, open_state = reduce_boundary(bd_active, k+sides, collisions, new_searchparams)
        end
        k += 1
    end
    # if error == 0
    #     p_length_id = min(Int64(ceil(total_t / d_length)), 32)
    #     N_frames = p_length_id*Int64(ceil(d_length/timestep))
    #     open(string(output_dir, "/totallength_", run, "_", p_length_id, "_animation.txt"), "a") do io
    #         for i in eachindex(positions_t)
    #             println(io, positions_t[i][1], " ", positions_t[i][2])
    #         end
    #         if length(positions_t) < N_frames
    #             for i in length(positions_t)+1:N_frames
    #                 println(io, ending[1], " ", ending[2])
    #             end
    #         end
    #     end
    # end
    # if error == 0
    #     open(string(output_dir, "/time_interval_", sides, "_", run, ".txt"), "a") do io
    #         for dist_hit in dt
    #             println(io, dist_hit)
    #         end
    #     end
    # end

    open(string(output_dir, "/final_pos_", sides, "_", run, ".txt"), "a") do io
        println(io, p.pos[1], " ", p.pos[2], " ", error, " ", k, " ", total_t)
    end

    return p, error, k
end



function non_selfavoiding_walk(bd::Billiard, reps::Int64)
    x0 = 3/4+1/16
    y0 = 3/16
    ϕ0 = 3*π/4
    
    for k in 1:20
        p = Particle(x0+k/200, y0+k/200, ϕ0)
        open(@sprintf("data/path_sinai_%d.txt",k),"w") do io
            println(io,p.pos[1], " ", p.pos[2])
        end
        for j in 1:reps
            i, t, ending, vel=bounce!(p,bd)
            open(@sprintf("data/path_sinai_%d.txt",k),"a") do io
                println(io,ending[1], " ", ending[2])
            end
        end
    end
end

function run_sim(reps::Int64, run::Int64)
    # savepaths = (string(output_dir,"/path_taken_", poly_name, "_", run, ".txt"), string(output_dir, "/final_pos_", poly_name, "_", run, ".txt"))
    # "data/path_taken_$poly_name_$run.txt"
    println("Run $run:")
    open(@sprintf("%s/starting_pos_%d_%d.txt",output_dir,sides,run), "w")
    # open(savepaths[2], "w")
    # open(savepaths[3], "w")
    # open("data/final_pos_$poly_name_$run.txt", "w")
    open(string(output_dir, "/final_pos_", sides, "_", run, ".txt"), "w")
    # open(string(output_dir, "/time_interval_", sides, "_", run, ".txt"), "w")
    for j in 1:reps
        bd = billiard_polygon(sides,1,[0,0])
        # bd_small = billiard_vertices([[0.,0.], (bd[1].ep + bd[1].sp)/2., bd[1].ep])
        x0, y0, ϕ0 = randominside_xyφ(bd)
        p = Particle(x0, y0, ϕ0)
        open(@sprintf("%s/starting_pos_%d_%d.txt",output_dir,sides,run), "a") do io
            println(io, x0, " ", y0, " ", ϕ0)
        end
        end_pos(bd,p)
        # open(savepaths[2], "a") do io
        #     println(io,p.pos[1], " ", p.pos[2], " ", error, " ", iterations)
        # end
        if j%10000 == 0
            println(j)
        end
    end
end

function edge_experiment(reps::Int64)
    savepaths = ("data/pathlength_triangle_edge.txt", "data/inflection_angles_triangle_edge.txt")
    open("data/final_pos_triangle_edge.txt","w")
    # open(savepaths[1], "w")
    # open(savepaths[2], "w")
    for j in 1:reps
        bd = billiard_polygon(sides, 1, [0,0])
        p = randominside(bd)
        p.vel = -1*p.vel
        i, t, cp = next_collision(p, bd)
        p.pos = cp
        p.vel = -1 * p.vel
        p, error, iterations = end_pos(bd, p, savepaths)
        if j %1000 == 0
            println(j)
        end
    end
end

# dict_polynames = Dict([(3, "triangle"), (4, "square"), (5, "pentagon"), (6, "hexagon"), (7, "septagon"), (8, "octagon"), (9, "nonagon"), (10, "decagon")])
global reps = parse(Int64, ARGS[1])
global sides = parse(Int64, ARGS[2])
global run = parse(Int64, ARGS[3])
global output_dir = ARGS[4]
# global poly_name = dict_polynames[sides]

# sides = 3    

# x0 = 0.24969026091209068
# y0 = -0.166066331105874
# ϕ0 = 2.67256614

# p1 = Particle(x0, y0, ϕ0)
# p2 = Particle(x0 + 0.01*cos(ϕ0), y0 + 0.01*sin(ϕ0), ϕ0)

# find_pretty(bd, p1, ("data/path1.txt", "data/positions1.txt"))
# find_pretty(bd, p2, ("data/path2.txt", "data/positions2.txt"))



run_sim(reps, run)

# non_selfavoiding_walk(billiard_polygon(sides, 1, [0,0]), reps)