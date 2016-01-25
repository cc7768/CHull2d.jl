#
# Convex Hull
#
immutable ConvexHull{T<:Real}
    extremepoints::Vector{Point{2, T}}
end

function intersect(ch1::ConvexHull, ch2::ConvexHull)
    error("Intersection between convex sets is not currently implemented")
end

#
# Line Segment
#
immutable LineSegment{T}
    # Points
    p1::Point{2, T}
    p2::Point{2, T}

    # Slope
    slope::Float64
end

function LineSegment(p1, p2)
    x1 = p1[1]
    y1 = p1[2]
    x2 = p2[1]
    y2 = p2[2]

    slope = (y2 - y1)/(x2 - x1)

    return Line(p1, p2, slope)
end

function evaluatey(x::T, l::LineSegment)
    x1 = l.p1[1]
    x2 = l.p2[1]

    if x1 < x < x2
        y = l.p1[2] + l.slope*(x-x1)
    else
        error("x is not on line segment")
    end

    return y
end

function evaluatex(y::T, l::LineSegment)
    y1 = l.p1[2]
    y2 = l.p2[2]

    if y1 < y < y2
        x = l.p1[1] + (y-y1)/l.slope
    else
        error("y is not on line segment")
    end

    return x
end

function intersect(l1::LineSegment, l2::LineSegment)
    if l1.slope == l2.slope
        warn("Same slope -- Do not intersect")
        status = -1
        intersection = Point(Inf, Inf)
        out = (status, intersection)
    end
    # Pull out points
    x1 = l1.p1[1]
    y1 = l1.p1[2]
    x2 = l2.p1[1]
    y2 = l2.p1[2]
    x3 = l1.p1[1]
    y3 = l1.p1[2]
    x4 = l2.p2[1]
    y4 = l2.p2[2]

    # Get intersection
    denom = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
    x_intersect = (x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*y4 - y3*x4)/denom
    y_intersect = (x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*y4 - y3*x4)/denom
    intersection = Point(x_intersect, y_intersect)

    # Check if intersection is on both intervals
    if (x1 < x_intersect < x2) && (x3 < x_intersect < x4)
        out = (1, intersection)
    # If not, return intersection but negative status
    else
        warn("Does not intersect on segment")
        out = (-2, intersection)
    end

    return out
end

immutable Quadrant{T}
    # The origin of quadrant (x, y)
    origin::Point{2, T}

    # The end point in x and y directions
    _end_x::T
    _end_y::T

    # The lines that make up quadrant
    _line_x = LineSegment{T}
    _line_y = LineSegment{T}
end

function Quadrant(origin::Point, _end_x, _end_y)
    # Get origin points
    x = origin[1]
    y = origin[2]

    # Create lines
    _line_x = Line(origin, Point(_end_x, y))
    _line_y = Line(origin, Point(x, _end_y))

    return Quadrant(origin, _end_x, _end_y, _line_x, _line_y)
end

function intersect(q1::Quadrant, q2::Quadrant)
    # Pull out points
    o1_x = q1.origin[1]
    o1_y = q1.origin[2]
    o2_x = q2.origin[1]
    o2_y = q2.origin[2]

    return nothing
end
