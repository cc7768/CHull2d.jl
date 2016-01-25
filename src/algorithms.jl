#=
Directory
---------
Akl Toussaint = line 13
Graham Scan = line 89
Monotone Chain = line 110
Quickhull
=#

#=
Akl Toussaint algorithm
=#
function _akltoussaint{T}(p::Vector{Point{2, T}})

    # Get number of points
    npts = length(p)
    if npts <= 4
        return p
    end

    # Get min and max values for x and y
    xmin, xminind, xmax, xmaxind = extremax(p)
    ymin, yminind, ymax, ymaxind = extremay(p)
    uptsind = unique([xminind, ymaxind, xmaxind, yminind])
    upts = p[uptsind]
    nu = length(uptsind)

    if nu == 4
        out = _akltoussaint4(p, upts, npts)
    elseif nu == 3
        out =  _akltoussaint3(p, upts, npts)
    else
        warn("Could not apply Akl Toussaint. Returning all points")
        out = p
    end

    return out
end

function _akltoussaint4{T<:Real}(p::Vector{Point{2, T}}, upts::Vector{Point{2, T}}, npts::Int64)

    # Get the points from upts
    pt1, pt2, pt3, pt4 = upts

    # Allocate space for output of each triangle
    out = Array(Bool, npts)

    # Create a "T inverse" matrix for each triangle
    tinv_1 = create_t_inv(pt1, pt2, pt3)
    tinv_2 = create_t_inv(pt1, pt4, pt3)

    @inbounds for n=1:npts
        temp = (p[n] - pt3)

        λ1_1, λ2_1 = _unsafe_mult_2by2_point(tinv_1, temp)
        λ1_2, λ2_2 = _unsafe_mult_2by2_point(tinv_2, temp)

        out[n] = check_λ(λ1_1, λ2_1) * check_λ(λ1_2, λ2_2)
    end

    return p[out]
end

function _akltoussaint3{T<:Real}(p::Vector{Point{2, T}}, upts::Vector{Point{2, T}}, npts::Int64)

    # Get the points from upts
    pt1, pt2, pt3 = upts

    # Allocate space for output of triangle
    out = Array(Bool, npts)

    # Create a "T inverse" matrix for each triangle
    tinv = create_t_inv(pt1, pt2, pt3)

    @inbounds for n=1:npts
        temp = (p[n] - pt3)

        λ1, λ2 = _unsafe_mult_2by2_point(tinv, temp)

        out[n] = check_λ(λ1, λ2)
    end

    return p[out]
end

#
# Graham Scan Algorithm
#
function _grahamscan{T<:Real}(p::Vector{Point{2, T}})

    # Get yminind
    xminind = indmin(p)
    xmin = p[xminind]

    # Create function to sort by angles
    lt(a, b) = xmin == a ? true : xmin == b ? false : ccw(xmin, a, b)
    psort = sort!(p, lt=lt)

    # Add the starting point at end so we go full circle
    push!(psort, xmin)

    # Call function that works around points
    ep = wrappoints(psort)

    return ep[1:end-1]
end

#
# Monotone Chain
#
function _monotonechain{T<:Real}(p::Vector{Point{2, T}})

    # Sort points
    psort = sort!(p)

    # Split into upper and lower
    xl, xu = split_xL_xU(psort)

    lh = wrappoints(xl)
    uh = wrappoints(reverse!(xu))

    return [lh[1:end-1]; uh[1:end-1]]
end

# function _quickhull{T<:Real}(p::Vector{Point{2, T}})

#     # Note: I suspect this algorithm can be implemented recursively
#     #       by calling it on two points, separating them into two
#     #       sets and then calling it on the endpoints of those
#     #       sets. Think carefully about how to implement this the
#     #       right way.
# end
