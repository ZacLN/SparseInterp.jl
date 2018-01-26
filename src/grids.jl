abstract type BasisFunction end
abstract type Linear  <: BasisFunction end
abstract type Quadratic <: BasisFunction end

struct NodeLocation
    l::UInt8
    j::UInt16
    dj::UInt16
end
NodeLocation(x::Float64) = NodeLocation(position(x)...)

mutable struct AdaptiveGrid
    nodeinfo
    active::BitArray{1}
    overlap::Vector{Vector{Bool}}
end

mutable struct NGrid{D, B}
    L::Vector{Int}
    bounds::Array{Float64, 2}
    grid::Array{Float64, 2}
    covers::Array{UInt16, 2}
    covers_dM::Array{UInt16, 2}
    covers_loc::Vector{UInt32}

    adapt::AdaptiveGrid

    IDs::Vector{Vector{Int}}
    Bs::Vector{Vector{Float64}}
end


# Returns the tensor grid of a vector of univariate grids
function Base.kron(X::Vector{Vector{Float64}})
    D = length(X)
    lxs = Int[length(x) for x in X]
    G = zeros(prod(lxs), D)
    s = 1
    for d = 1:D
        snext = s * lxs[d]
        for j = 1:prod(lxs)
            G[j, d] = X[d][div(rem(j - 1, snext), s) + 1]
        end
        s = snext
    end
    return G
end

Base.kron(x::Vector{Float64}) = x

"""
    TensorGrid(L)

Returns a tensor grid of univariate grids of levels given by L = [l1,l2,...]
"""
function TensorGrid(L::Vector{Int})
    G = ndgrid(Vector{Float64}[cc_g(i) for i in L]...)
    G = hcat([vec(g) for g in G]...)
end


"""
    SmolyakSize(L[, mL])

Returns the number of nodes in levels mL of a Smolyak grid of maximum
of maximum level L = [l1,l2...]
"""
function SmolyakSize(L::Vector{Int}, mL::UnitRange{Int} = 0:maximum(L))
    D = length(L)
    m = Int[dM(l) for l = 1:maximum(L) + D]
    S = 0
    for l = mL
        for covering in comb(D, D + l)
            if all(covering .≤ L + 1)
                s = m[covering[1]]
                for i = 2:length(covering)
                    s *= m[covering[i]]
                end
                S += s
            end
        end
    end
    return S
end


"""
    SmolyakGrid(L[, mL])

 juli
"""
function SmolyakGrid(L::Vector{Int}, mL::UnitRange{Int} = 0:maximum(L))
    D = length(L)
    dg = Vector{Float64}[cc_dg(l) for l = 1:maximum(L) + D]
    G = Array{Array{Float64}}(0)
    index = Array{Array{Int}}(0)
    for l = mL
        for covering in comb(D, D + l)
            if all(covering .≤ L + 1)
                push!(G, kron(dg[covering]))
                push!(index, repmat(covering', size(G[end], 1)))
            end
        end
    end
    # G = vcat(G...)::Array{Float64,2}
    # index = vcat(index...)::Array{Int,2}
    # return G,index
    return (vcat(G...), vcat(index...))::Tuple{Array{Float64, 2}, Array{Int, 2}}
end

"""
    level(x)

Computes the minimum level of a point
"""
function level(x::Float64)
    l = 0
    if x == 0.5
        l =  1
    elseif x == 0.0 || x == 1.0
        l = 2
    else
        for l = 3:12
            mod(x, 1 / 2^(l - 1)) == 0.0 && break
        end
    end
    return l
end
level(G::NGrid) = level(G.grid)
level(X::Array{Float64}) = vec(sum(map(level, X), 2) - size(X, 2))


"""
    NGrid(L,bounds,B)

# Construction
Construct a Smolyak grid.
The vector of integers L determines the maximum level of the grid in each
dimension. The bounds can be optionally changed from the default of [0,1]^D and
basis function B can be either Linear or Quadratic.

# Interpolation
Grid objects are callable taking two arguments. The first is a vector containing
function values at the grid nodes. The second array contains rows of points at
which the interpolant is to be evaluated.
"""
function NGrid(L::Vector{Int}, bounds::Array{Float64} = [0, 1] .* ones(1, length(L));B::Type{BT} = Linear) where {BT <: BasisFunction}
    grid, ind = SmolyakGrid(L)
	   covers = UInt16.(unique(ind, 1))
    covers_loc = zeros(UInt32, size(covers, 1))
    hind = vec(mapslices(hash, ind, 2))
    hcovers = vec(mapslices(hash, covers, 2))
    for i = 1:size(covers, 1)
        covers_loc[i] = findfirst(hind, hcovers[i])
    end

    G = NGrid{length(L), B}(L,
                        bounds,
                        grid,
                        covers,
                        dM.(covers),
                        covers_loc,
                        AdaptiveGrid([], zeros(Bool, size(grid, 1)), []),
                        Vector{Int}[Int[] for i = 1:size(grid, 1)],
                        Vector{Float64}[Float64[] for i = 1:size(grid, 1)]
                        )
    G.adapt.active = (level(G) .== maximum(level(G)))
    buildW(G, hind, hcovers)
    return G
end

Base.show(io::IO, G::NGrid) = println(io, typeof(G), ": $(size(G.grid,1))pts")
Base.length(G::NGrid) = size(G.grid, 1)
Base.size(G::NGrid) = size(G.grid)
Base.values(G::NGrid) = nUtoX(G.grid, G.bounds)


Base.getindex(G::NGrid, args...) = getindex(G.grid, args...)
