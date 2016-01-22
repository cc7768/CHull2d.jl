#
# General utils
#
@inline orientation(a::Point, b::Point, c::Point) =
    (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x)
cw(a, b, c) = orientation(a, b, c) < 0.0 ? true : false
ccw(a, b, c) = orientation(a, b, c) >= 0.0 ? true : false

#
# Utils for Akl Toussaint
#
"""
Creates the "T inverse" matrix described by wikipedia and
above that will right multiply points to get the Barycentric
coordinates
"""
function create_t_inv(a::Point, b::Point, c::Point)
    # Pull out points
    a1, a2 = a
    b1, b2 = b
    c1, c2 = c

    # Build t inverse matrix
    temp = -a2* b1 + a1 *b2 + a2 *c1 - b2 *c1 - a1 *c2 + b1 *c2
    return [(b2-c2)/temp (c1-b1)/temp; (c2-a1)/temp (a1-c1)/temp]
end

# One liner to check if λ values are in (0, 1) and sum to less than 1
@inline check_λ(λ1, λ2) = (0 < λ1 < 1) && (0 < λ2 < 1) && (λ1 + λ2 < 1) ? false : true

"""
Determines which points in a vector of points lie inside or outside the
triangle determined by the triangle abc.
"""
function pt_out_tri{T}(x::Vector{Point{T}}, a::Point, b::Point, c::Point)
    # Get number of points in x
    npts = length(x)

    # Allocate space
    out = Array(Bool, npts)

    # Create our "T inverse" matrixs
    tinv = create_t_inv(a, b, c)

    for n=1:npts
        # Get the vector we need to multiply by T inverse to get
        # the Barycentric coordinates
        temp = x[n] - c

        # Get Barycentric coordinates
        λ1, λ2 = tinv*temp
        out[n] = check_λ(λ1, λ2)
    end

    return out
end

#
# Wrap around set of points
#
"""
Takes a vector of points (which must be sorted either by
relative  or x-value) and wraps around them by assuring
that they only make counter-clockwise turns.
"""
function wrappoints{T<:Real}(p::Vector{Point{T}})

    # Get number of points
    npts = length(p)

    # Put first two points inside the convex hull
    ch = [p[1], p[2]]
    sizehint!(ch, npts)

    # Starting value for number of extreme points
    m = 2

    # Will check each point
    @inbounds for i=3:npts
        # If it makes a clockwise turn then the
        # last element of the convex hull doesn't
        # belong, so pop it out. Make sure that at
        # least two elements are always in it
        while m>1 && cw(ch[m-1], ch[m], p[i])
            pop!(ch)
            m = m-1
        end
        push!(ch, p[i])
        m = m+1
    end

    return ch
end

#
# Utils for Monotone Chain
#
function split_xL_xU{T<:Real}(psort::Vector{Point{T}})
    # Allocate space to identify where points belong
    # We make endpoints true because we want them in both sets
    npts = length(psort)
    xL_bool = Array(Bool, npts)
    xL_bool[1] = true; xL_bool[end] = true
    xU_bool = Array(Bool, npts)
    xU_bool[1] = true; xU_bool[end] = true

    # Get the endpoints and slope
    pmin = psort[1]; pmax = psort[end]
    xmin = pmin.x; ymin = pmin.y
    slope = (pmax.y - pmin.y)/(pmax.x - pmin.x)

    # Iterate through all points except endpoints
    @inbounds for i=2:npts-1
        # Pull out x and y values for ith point
        p_i = psort[i]
        xi, yi = p_i.x, p_i.y

        # Compute where y on line would be
        yi_line = ymin +  slope*(xi - xmin)

        if yi_line <= yi
            xL_bool[i] = false
            xU_bool[i] = true
        else
            xL_bool[i] = true
            xU_bool[i] = false
        end
    end
    xL = psort[xL_bool]
    xU = psort[xU_bool]

    return xL, xU
end

