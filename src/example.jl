using FixedSizeArrays
using PlotlyJS

import Base: *, isless

include("utils.jl");
include("algorithms.jl");

p = [Point(randn(), randn()) for i=1:250];
ep = _akltoussaint(p)

# a = Point(0.0, 0.0)
# b = Point(1.0, 0.1)
# c = Point(0.5, 0.4)
# d = Point(1.0, 1.0)
# e = Point(0.1, 1.0)

# wrappoints([a, b, c, d, e, a])

ep_g = _grahamscan(p)
ep_m = _monotonechain(p)

# Plot all points
t1 = scatter(;x=[i[1] for i in p], y=[i[2] for i in p], mode="markers")

# Plot Convex hull from graham scan
t2 = scatter(;x=[i[1] for i in ep_g], y=[i[2] for i in ep_g], mode="lines")

# Plot upper hull
t3 = scatter(;x=[i[1] for i in ep_m], y=[i[2] for i in ep_m], mode="lines")

show(Plot([t1, t2, t3]))
