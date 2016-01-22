#=
This file contains a few algorithms that are used to compute convex
hulls in two dimensions.

Currently included are:

* Graham Scan Algorithm
* Monotone Chain Algorithm

Algorithms that should be added are:

* Quickhull Algorithm
* Chan Algorithm
* ?
=#

# Load dependencies
using FixedSizeArrays

# Base methods
import Base: *, isless, isequal

# Create a type to hold convex hulls
immutable ConvexHull{T<:Real}
    points::Vector{Point{2, T}}
    extremepoints::Vector{Point{2, T}}
end

#
# Include all algorithms
#
include("algorithms.jl")
include("utils.jl")

const ALGDICT = Dict{Symbol, Function}()
ALGDICT[:MonotoneChain] = _monotonechain
ALGDICT[:GrahamScan] = _grahamscan
# ALGDICT[:QuickHull] = _quickhull

"""
Creates the convex hull of a set of points

* points : The set of points that you would like to take the convex hull of
* algorithm : Which convex hull algorithm to use. Currently implemented are
              [:GrahamScan, :MonotoneChain]
* _at : Whether or not to use the Akl Toussaint pruning
"""
function ConvexHull(points::Vector{Point}, algorithm=:MonotoneChain, _at=true)
    # First prune points
    p = _at ? _akltoussaint(points) : points

    # Apply algorithm
    ep = ALGDICT[algorithm](p)

    return ConvexHull(points, ep)
end
