using Benchmark
using FixedSizeArrays

include("utils.jl");
include("algorithms.jl");

function graham_scan_vanilla(p::Vector{Point{2, Float64}})
    ep = _grahamscan(p)

    return ep
end

function graham_scan_wat(p::Vector{Point{2, Float64}})
    p = _akltoussaint(p)
    ep = _grahamscan(p)

    return ep
end

function monotonechain_vanilla(p::Vector{Point{2, Float64}})
    ep = _monotonechain(p)

    return ep
end

function monotonechain_wat(p::Vector{Point{2, Float64}})
    p = _akltoussaint(p)
    ep = _monotonechain(p)

    return ep
end

srand(1112016)
npoints = [25, 250, 2500, 25000, 250000, 2500000]
nreps = [1000, 500, 100, 50, 25, 5]

for i=1:length(npoints)

    # Get parameters
    np = npoints[i]
    nr = nreps[i]

    p = [Point(randn(), randn()) for j=1:np]

    gsv() = graham_scan_vanilla(p)
    gsat() = graham_scan_wat(p)
    mcv() = monotonechain_vanilla(p)
    mcat() = monotonechain_wat(p)

    println("With $np points")
    println(compare([gsv, gsat, mcv, mcat], nr))
end

