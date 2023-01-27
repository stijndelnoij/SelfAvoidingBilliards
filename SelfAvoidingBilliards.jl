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
include("timeevolution/boundary_reduction.jl")


end#module
