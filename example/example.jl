using FixedSizeArrays
using PlotlyJS
include("../src/CHull2D.jl")

#
# LineSegment Intersection
#

# Simple example
p1 = Point(0.0, 0.0)
p2 = Point(0.5, 0.5)
p3 = Point(0.0, 0.5)
p4 = Point(0.5, 0.0)
l1 = CHull2D.LineSegment(p1, p2)
l2 = CHull2D.LineSegment(p3, p4)
println("Should be true", ":", intersect(l1, l2))

# TODO: Should be able to handle Inf instead of 1e12
p5 = Point(0.0, 0.25)
p6 = Point(1e12, 0.25)
l3 = CHull2D.LineSegment(p5, p6)
println("Should be true", ":", intersect(l1, l3))

#
# Convex Hull with Quadrant
#
p = [Point(randn(), randn()) for i=1:25]
ch = CHull2D.convexhull(p)
q = CHull2D.Quadrant(Point(0.0, 0.0), Inf, Inf)
chiq = intersect(ch, q)

t1 = scatter(;x=[foo[1] for foo in chiq], y=[foo[2] for foo in chiq])
t2 = scatter(;x=[foo[1] for foo in ch.extremepoints], y=[foo[2] for foo in ch.extremepoints])
Plot([t1, t2])
