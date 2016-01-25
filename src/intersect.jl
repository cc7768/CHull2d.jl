# Intersection methods
# TODO: Change intersection methods -- This breaks with infinite slope
function intersect(l1::LineSegment, l2::LineSegment)

    if l1.slope == l2.slope
        warn("Same slope -- Do not intersect")
        intersection = Point(Inf, Inf)
        out = (false, intersection)
    end

    # Pull out points
    x1 = l1.p1[1]
    y1 = l1.p1[2]
    x2 = l1.p2[1]
    y2 = l1.p2[2]
    x3 = l2.p1[1]
    y3 = l2.p1[2]
    x4 = l2.p2[1]
    y4 = l2.p2[2]

    # Get intersection
    denom = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
    x_intersect = (x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*y4 - y3*x4)/denom
    y_intersect = (x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*y4 - y3*x4)/denom
    intersection = Point(x_intersect, y_intersect)

    # Check if intersection is on both intervals
    if in(x_intersect, Interval(x1, x2)) && in(x_intersect, Interval(x3, x4))
        out = (true, intersection)
    # If not, return intersection but negative status
    else
        out = (false, intersection)
    end
in
    return out
end

function intersect(ch1::ConvexHull, ch2::ConvexHull)
    error("Intersection between convex sets is not currently implemented")
end

function intersect(ch::ConvexHull, q::Quadrant)

    # Get information about the set of extreme points
    ext_pts = ch.extremepoints
    xmin, _, xmax, _ = extremax(ext_pts)
    ymin, _, ymax, _ = extremay(ext_pts)
    npts = length(ch)

    # Get information about quadrant
    o = q.origin
    lx = q.lx
    ly = q.ly

    # Start a convex hull to store outpoints in
    prev_pt = ext_pts[end]
    outhull = in(prev_pt, q) ? [prev_pt] : []

    # Will check between each point individually
    # Points need to be ordered either clockwise
    # or counter-clockwise
    for i=1:npts
        # Current extreme point and corresponding line
        c_ep = ext_pts[i]
        c_ls = LineSegment(prev_pt, c_ep)

        # Compute intersections with lines
        stat_x, c_inter_x = intersect(lx, c_ls)
        stat_y, c_inter_y = intersect(ly, c_ls)

        # If there is intersection, push to hull
        stat_x ? push!(outhull, c_inter_x) : nothing
        stat_y ? push!(outhull, c_inter_y) : nothing

        # Add current point if it is in quadrant
        in(c_ep, q) ? push!(outhull, c_ep) : nothing
        prev_pt = c_ep

    end

    return outhull
end


# function intersect(ch::ConvexHull, q::Quadrant)
#     # TODO: Need to consider special cases where it is not binding on one side
#     # but is on the other. It should be as simple as returning all extrema to
#     # the right side of the binding extreme. Right thing to do here is break
#     # this function into cases and write the cases out, that will make it
#     # look much less messy

#     # Get information about the set of extreme points
#     ext_pts = ch.extremepoints
#     _, xmin, _, xmax = extremax(ext_pts)
#     _, ymin, _, ymax = extremay(ext_pts)
#     npts = length(ep)

#     # Get information about quadrant
#     o = q.origin

#     # Counters to determine whether points belong in the
#     # set or not. vcount and hcount are for number of intersections
#     # the line makes with the convex hull each way. If either are
#     # even then origin is not in intersection
#     checknext = true
#     storecurr = false
#     vcount = 0
#     hcount = 0
#     # Start a convex hull to store outpoints in
#     outhull = []

#     # Will check between each point individually
#     # Points need to be ordered either clockwise
#     # or counter-clockwise
#     for i=2:npts+1
#         # Pull out points
#         if i<=npts
#             x1 = ep[i-1]
#             x2 = ep[i]
#         else
#             x1 = ep[npts]
#             x2 = ep[1]
#         end

#         # First add point if we are supposed to store
#         if storecurr
#             push!(outhull, x1)
#         end

#         # Check whether intersects
#         w1intersect = evalliney(x1, x2, o.x)
#         w2intersect = evallinex(x1, x2, o.y)
#         cond1 = ((x2.x - o.x)*(x1.x - o.x) < 0) && (w1intersect >= o.y)
#         cond2 = ((x2.y - o.y)*(x1.y - o.y) < 0) && (w2intersect >= o.x)

#         # Check first condition
#         if cond1
#             storecurr = !storecurr
#             push!(outhull, Point(w1, w1intersect))
#             vcount += 1
#         # Account for if point is right on top of value
#         elseif isapprox(x2.x, o.x) && isapprox(x2.y, w1intersect)
#             if storecurr
#                 push!(outhull, x2)
#             end
#             storecurr = !storecurr
#             vcount += 1
#         end

#         # Check second condition
#         if cond2
#             storecurr = !storecurr
#             # Make sure we don't double count if both conditions hold
#             push!(outhull, Point(w2intersect, w2))
#             hcount += 1
#         elseif isapprox(x2.y, o.y) && isapprox(x2.x, w2intersect)
#             if storecurr
#                 push!(outhull, x2)
#             end
#             storecurr = !storecurr
#             hcount += 1
#         end
#     end

#     if (hcount*vcount)%2!=0 || (o.x<minx(ep)[1].x && o.y<miny(ep)[1].y)
#         push!(outhull, o)
#     end

#     return outhull
# end
