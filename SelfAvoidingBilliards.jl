__precompile__()

"""
A Julia package for self avoiding billiard systems in two dimensions.

It provides a flexible, easy-to-use, and intuitive framework for
fast implementation of billiards of arbitrary construction.

This code is (almost) a complete copy of the dynamical billiards Julia package, with some adjustments to accomodate self avoidance.
"""
module SelfAvoidingBilliards

using LinearAlgebra
using StaticArrays
using Elliptic
using DataFrames
import Base: show, eltype, getindex

const SV = SVector{2}
export SVector

cossin(a) = ((x, y) = sincos(a); (y, x))

##########################################
# Core                                   #
##########################################
include("billiards/particles.jl")
include("billiards/obstacles.jl")
include("billiards/billiardtable.jl")
include("billiards/standard_billiards.jl")

include("timeevolution/collisions.jl")
include("timeevolution/propagation.jl")
include("timeevolution/timeseries.jl")
include("timeevolution/highleveltimes.jl")



###################
# Update messages #
###################
using Scratch
display_update = true
version_number = "4"
update_name = "update_v$(version_number)"

function __init__()
if display_update
    # Get scratch space for this package
    versions_dir = @get_scratch!("versions")
    if !isfile(joinpath(versions_dir, update_name))
        printstyled(
            stdout,
            """
            \nUpdate message: DynamicalBilliards v$(version_number)
            Plotting & animating has moved entirely to the Makie ecosystem.
            It is now provided by InteractiveDynamics.jl.
            See the documentation online for the new API.
            """;
            color = :light_magenta,
        )
        touch(joinpath(versions_dir, update_name))
    end
end
end


end#module
