#
# Convex Hull
#
immutable ConvexHull{T}
    extremepoints::Vector{Point{2, T}}

    # TODO: Add inner constructor
end

length(ch::ConvexHull) = length(ch.extremepoints)
isequal(ch1::ConvexHull, ch2::ConvexHull) = ch1.extremepoints == ch2.extremepoints

#
# Interval
#
immutable Interval
    lb
    ub

    function Interval(lb, ub)
        if lb<ub
            return new(lb, ub)
        else
            return new(ub, lb)
        end
    end
end

in(x, i::Interval) = (i.lb <= x <= i.ub) ? true : false

#
# Line Segment
#
immutable LineSegment{T}
    # Points
    p1::Point{2, T}
    p2::Point{2, T}

    # Slope
    _slope::Float64
    _norm::Float64
end

function LineSegment(p1, p2)
    x1 = p1[1]
    y1 = p1[2]
    x2 = p2[1]
    y2 = p2[2]

    _slope = (y2 - y1)/(x2 - x1)
    _norm = sqrt((y2-y1)*(y2-y1) + (x2-x1)*(x2-x1))

    return LineSegment(p1, p2, _slope, _norm)
end

isfinite(l::LineSegment) = all((isfinite(l.p1[1]), isfinite(l.p1[2]), isfinite(l.p2[1]), isfinite(l.p2[2])))

function evaluatey(l::LineSegment, x::Float64)
    x1, y1 = l.p1
    x2, y2 = l.p2

    if x1 <= x <= x2
        y = y1 + l._slope*(x-x1)
    else
        error("x is not on line segment")
    end

    return y
end

function evaluatex(l::LineSegment, y::Float64)
    y1 = l.p1[2]
    y2 = l.p2[2]

    if y1 < y < y2
        x = l.p1[1] + (y-y1)/l._slope
    else
        error("y is not on line segment")
    end

    return x
end

# TODO: This isn't quite right. It returns distance 0 for points on same line
# (but not in line segment)
function distance(l::LineSegment, p::Point)
    # Pull out points
    x0, y0 = p[1], p[2]
    x1, y1 = l.p1[1], l.p1[2]
    x2, y2 = l.p2[1], l.p2[2]

    # Compute numerator and denominator of Wikipedia expression
    numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
    denominator = l._norm

    return numerator/denominator
end


#
# Quadrant
#
immutable Quadrant{T}
    # The origin of quadrant (x, y)
    origin::Point{2, T}

    # The end point in x and y directions (should be Inf or -Inf)
    lx::LineSegment
    ly::LineSegment

    # Functions which say whether points are in correct piece of quadrant
    dir_x::Function
    dir_y::Function

end

function Quadrant(origin::Point, _end_x, _end_y)
    # Get origin points
    x = origin[1]
    y = origin[2]

    # If end is positive infinity then want to check right of origin
    # TODO: Come back and allow this to be Inf -- The intersection
    # method for line segments can't handle this currently
    if _end_x == Inf
        dir_x(p) = p[1] >= x
        lx = LineSegment(origin, Point(1e12, y))
    elseif _end_x == -Inf
        dir_x(p) = p[1] <= x
        lx = LineSegment(Point(-1e12, y), origin)
    end

    # If end is positive infinity then want to check above origin
    if _end_y == Inf
        dir_y(p) = p[2] >= y
        ly = LineSegment(origin, Point(x, Inf))
    elseif _end_y == -Inf
        dir_y(p) = p[2] <= y
        ly = LineSegment(Point(x, -Inf), origin)
    end

    return Quadrant(origin, lx, ly, dir_x, dir_y)
end

in(x::Point, q::Quadrant) = (q.dir_x(x) && q.dir_y(x)) ? true : false

