#=
This file contains a few algorithms that are used to compute convex
hulls in two dimensions.

Currently included are:

* Graham Scan Algorithm
* Monotone Chain Algorithm
* Quickhull Algorithm

Algorithms that should be added are:

* Chan
* ?
=#

# It turns out to be useful to have our own
# point type so that we can iterate in 1
# dimension over points.
include("point.jl")


# Create a type to hold convex hulls
immutable ConvexHull{T<:Real}
    points::Vector{Point{T}}
    extremepoints::Vector{Point{T}}
end

#
# Include all algorithms
#
include("algorithms.jl")
include("utils.jl")

const ALGDICT = Dict{Symbol, Function}()
ALGDICT[:MonotoneChain] = _monotonechain
ALGDICT[:GrahamScan] = _grahamscan
ALGDICT[:QuickHull] = _quickhull

function ConvexHull(points::Vector{Point}, algorithm=:MonotoneChain)
    # First prune points
    p = _akltoussaint(points)

    # Apply algorithm
    ep = ALGDICT[algorithm](p)

    return ConvexHull(points, ep)
end
