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
module CHull2D

# Load dependencies
using FixedSizeArrays

# Base methods
import Base: *, in, intersect, isfinite, isless, isequal, length

# Include file of all types
include("CHullTypes.jl")
include("utils.jl")

#
# Include all algorithms
#
include("intersect.jl")
include("algorithms.jl")

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
function convexhull{T}(points::Vector{Point{2, T}}; algorithm::Symbol=:MonotoneChain, _at::Bool=true)
    # TODO: Should make sure that there are at least 3 not co-linear points
    # First prune points
    points = _at ? _akltoussaint(points) : points

    # Apply algorithm
    ep = ALGDICT[algorithm](points)

    return ConvexHull(ep)
end

export ConvexHull, LineSegment, Quadrant, Point, # Types
       convexhull, wrappoints, pt_out_tri,
       evaluatex, evaluatey, ccw, cw, orientation

end

