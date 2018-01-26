function addcover!(G::NGrid{D, B}, L::Vector{Int}) where {D,B}
    any(all(G.covers .== L', 2)) && return
    bf = B == Linear ? cc_bf_l : cc_bf_q
    x   = kron([SparseInterp.cc_dg(i) for i in L])
    ind = repmat(map(Int16, L'), size(x, 1))
    level_M = map(i -> Int16(SparseInterp.M(level(i))), G.grid)

    for id = 1:size(x, 1)
        push!(G.Bs, [])
        push!(G.IDs, Int[])
        for i = 1:length(G)
            b = 1.0
            for d = 1:length(G.L)
                b *= bf(x[id, d], G.grid[i, d], level_M[i, d])
            end
            b > 0 && push!(G.Bs[end], -b)
            b > 0 && push!(G.IDs[end], i)
        end
        push!(G.Bs[end], 1.0)
        push!(G.IDs[end], length(G) + id)
    end
    G.grid = [G.grid; x]
    ind = map(level, G.grid)

    G.covers = [G.covers;L']
    G.covers_loc = Int32[findfirst(all(G.covers[i:i, :] .== ind, 2)) for i = 1:size(G.covers, 1)]
    G.covers_dM = map(x -> SparseInterp.dM(Int(x)), G.covers)
    G.adapt.active = [G.adapt.active;zeros(Bool, size(x, 1))]
    G.L = vec(maximum(G.covers, 1)) - 1
    return
end

function Base.sort!(G::NGrid)
    lev = level(G)
    try
        for i = 2:length(G)
            @assert lev[i - 1] ≤ lev[i]
        end
    catch
        id = zeros(Int, 0)
        for l = 0:maximum(lev)
            id = [id;find(lev .== l)]
        end
        G.grid = G.grid[id, :]
        ind = map(level, G.grid)
        G.covers = map(UInt16, unique(ind, 1))
        G.covers_dM = map(x -> dM(Int(x)), G.covers)
        G.covers_loc = Int32[findfirst(all(G.covers[i:i, :] .== ind, 2)) for i = 1:size(G.covers, 1)]
        SparseInterp.buildW(G)
    end
end

function iscoverfull(G::NGrid, L::Vector{Int})
    @assert length(L) == length(G.L)
    n = 0
    D = length(L)
    for i = 1:length(G)
        b = true
        for j = 1:D
            b = b && level(G.grid[i, j]) == L[j]
        end
        b && (n += 1)
    end
    return n == prod(map(dM, L))
end

looseid(G::NGrid) = G.covers_loc[end] + prod(G.covers_dM[end, :]):size(G.grid, 1)



"""
    grow!(G,id,bounds)

Adapt grid G at node number 'id'. 2*D nodes are added (except when G[id] is on
an edge in which case 2*D-e nodes are added where e is the number of edges G[id]
sits on). 'bounds' limits the maximum level the grid can grow to.
"""
function grow!(G::NGrid{D, BF}, id::Int, bounds::Vector{Int} = 12 * ones(Int, length(G.L))) where {D,BF}
    !G.adapt.active[id] && (return nothing)
    G.adapt.active[id] = false

    x       = G.grid[id, :]
    n       = D * 2 - sum((x[d] == 0.0 || x[d] == 1.0) for d = 1:D)
    Lid     = sum(map(level, G.grid[id, :])) - D
    X, ind   = SmolyakGrid(min(map(level, x) + 1, bounds), Lid + (1:1))

    id1     = (hX = mapslices(hash, G.grid, 2);Bool[!in(hash(X[i, :]), hX) for i = 1:size(X, 1)])
    all(!id1) && (return nothing)
    X = X[id1, :]
    id1 = sortperm(Float64[norm(X[i, :] - vec(x)) for i = 1:size(X, 1)])[1:min(n, size(X, 1))]
    X = X[id1, :]
    ind = map(level, X)

    hind = mapslices(hash, ind, 2)
    firstloc = length(G)
    for i = 1:n
        Hind = mapslices(hash, map(level, G.grid), 2) # need not be completely recalculated, can insert new hash
        if in(hind[i], Hind)
            insertloc = findlast(Hind, hind[i])
        else
            insertloc = max(findlast(sum(map(level, G.grid), 2), sum(ind[i, :])), length(G))
        end
        firstloc = min(firstloc, insertloc)

        G.grid          = [G.grid[1:insertloc, :];X[i, :]'G.grid[insertloc + 1:end, :]]
        G.Bs            = [G.Bs[1:insertloc];[Float64[]];G.Bs[insertloc + 1:end]]
        G.IDs           = [G.IDs[1:insertloc];[Float64[]];G.IDs[insertloc + 1:end]]
        G.adapt.active  = [G.adapt.active[1:insertloc];true;G.adapt.active[insertloc + 1:end]]
    end

    G.adapt.overlap = getparents(G)
    rebuildW(G, firstloc:length(G))

    return
end

function grow!(G::NGrid, id::Vector{Int}, bounds::Vector{Int} = 12 * ones(Int, length(G.L)))
    X = [G.grid[i, :] for i in id]
    for x in X
        grow!(G, findrow(G.grid, x), bounds)
    end
end



function rebuildW(G::NGrid{D, BF}, idr) where {D,BF}
    bf = (BF == Linear ? cc_bf_l : cc_bf_q)
    level_M = map(i -> Int16(SparseInterp.M(level(i))), G.grid)
    for j = idr
        G.Bs[j] = Float64[]
        G.IDs[j] = Int[]
        for i = 1:j - 1
            b = 1.0
            for d = 1:D
                b *= bf(G.grid[j, d], G.grid[i, d], level_M[i, d])
            end
            if b > 0
                push!(G.Bs[j], -b)
                push!(G.IDs[j], i)
            end
        end
        push!(G.Bs[j], 1.0)
        push!(G.IDs[j], length(G) + j)
    end
end

mutable struct box{D}
    lower::Vector{Float64}
    upper::Vector{Float64}
end

function box(G::NGrid{D, BF}, i::Int) where {D,BF}
    lower = Array{Float64}(D)
    upper = Array{Float64}(D)
    for d = 1:D
        p = level(G.grid[i, d])
        lower[d] = clamp(G.grid[i, d] - 1 / 2^(p - 1), 0, 1) + eps(Float64)
        upper[d] = clamp(G.grid[i, d] + 1 / 2^(p - 1), 0, 1) - eps(Float64)
    end
    return box{D}(lower, upper)
end

function Base.intersect(a::box{D}, b::box{D}) where {D}
    nool = a.lower[1] > b.upper[1] ||
           b.lower[1] > a.upper[1]
    for d = 2:D
        nool = nool ||
           a.lower[d] > b.upper[d] ||
           b.lower[d] > a.upper[d]
    end
    return !nool
end

"""
   getparents(G)

For each node of grid G this generates a Boolean Vector indicating whether the
basis function overlaps with that of the other nodes.
"""
function getparents(G::NGrid)
    T = [zeros(Bool, length(G)) for i = 1:length(G)]
    bx = [box(G, i) for i = 1:length(G)]
    @threads for i = 1:length(G)
        for j = 1:length(G)
            T[i][j] = intersect(bx[i], bx[j])
        end
    end
    return T
end


# Generate interpolation functions for single function approximation
for b in [(Linear, Lbj, cc_bf_l), (Quadratic, Qbj, cc_bf_q)]
    for D = 2:15
        coverloop = :(@inbounds for ii in nc
            b   = B[G.covers[ii, $D], $D] * B[G.covers[ii, 1], 1]
            id1 = J[G.covers[ii, $D], $D]
        end)
        for d in D - 1:-1:2
            push!(coverloop.args[2].args[2].args, :(b  *= B[G.covers[ii, $d], $d]))
            push!(coverloop.args[2].args[2].args, :(id1 = id1 * G.covers_dM[ii, $d] + J[G.covers[ii, $d], $d]))
        end

        push!(coverloop.args[2].args[2].args, :(id1 = id1 * G.covers_dM[ii, 1] + G.covers_loc[ii] + J[G.covers[ii, 1], 1];yi += b * w[id1]))

        f = :(function jl_ainterp(G::NGrid{$D, $(b[1])}, A::Vector{Float64}, xi::Array{Float64, 2}, y = zeros(Float64, size(xi, 1)))
            w       = getW(G, A)
            x       = nXtoU(xi, G.bounds)
            nx      = size(x, 1)
            nc      = 1:size(G.covers, 1)
            nG      = length(G)
            ns      = (G.covers_loc[end] + prod(G.covers_dM[end, :]))
            nc2     = ns:UInt32(nG)
            dr1 = 1:$D
            mL      = maximum(G.L) + 1
            J         = zeros(Int, mL, $D)
            B         = ones(mL, $D)
            # dM = map(x -> M(level(x)), G.grid)

            @threadsfixed [J, B] for i = 1:nx
                $(b[2])
                yi = 0.0
                coverloop

                i1, i2 = UInt32(1), UInt32(1)
                iswitch = false

                @inbounds for ii = nc2
                    if G.adapt.overlap[i1][ii] && G.adapt.overlap[i2][ii]
                        b = 1.0
                        for d = 1:$D
                            # b *= $(b[3])(x[i, d], G.grid[ii, d], dM[ii, d])
                            b *= $(b[3])(x[i, d], G.grid[ii, d], M(level(G.grid[ii, d])))
                            b == 0 && break
                        end
                        if b > 0
                            (yi += w[ii] * b)
                            iswitch = !iswitch
                            iswitch ? (i1 = ii::UInt32) : (i2 = ii::UInt32)
                        end
                    end
                end
                y[i] = yi
            end
            y
        end)

        subs!(f, :coverloop => coverloop)
        eval(current_module(), f)
    end
end

# Generate interpolation functions for multiple function approximation
for b in [(Linear, Lbj, cc_bf_l), (Quadratic, Qbj, cc_bf_q)]
    for D = 2:15
        coverloop = quote
            fill!(yi, 0.0)
            @inbounds for ii in nc
                b   = B[G.covers[ii, $D], $D] * B[G.covers[ii, 1], 1]
                id1 = J[G.covers[ii, $D], $D]
            end
        end
        for d in D - 1:-1:2
            push!(coverloop.args[4].args[2].args[2].args, :(b  *= B[G.covers[ii, $d], $d]))
            push!(coverloop.args[4].args[2].args[2].args, :(id1 = id1 * G.covers_dM[ii, $d] + J[G.covers[ii, $d], $d]))
        end
        push!(coverloop.args[4].args[2].args[2].args, :(id1 = id1 * G.covers_dM[ii, 1] + G.covers_loc[ii] + J[G.covers[ii, 1], 1]))


        push!(coverloop.args[4].args[2].args[2].args, :(for d = 1:nA yi[d] += b * w[id1, d] end))
        push!(coverloop.args, :(@inbounds for d = 1:nA y[i, d] = yi[d] end))


        f = :(function jl_ainterp(G::NGrid{$D, $(b[1])}, A::Array{Float64, 2}, xi::Array{Float64, 2}, y = zeros(Float64, size(xi, 1), size(A, 2)))
            w       = getW(G, A)
            x       = nXtoU(xi, G.bounds)
            nx      = size(x, 1)
            nA      = size(A, 2)
            nc      = 1:size(G.covers, 1)
            nG      = length(G)
            ns      = (G.covers_loc[end] + prod(G.covers_dM[end, :]))
            nc2     = ns:UInt32(nG)
            dr1 = 1:$D
            mL      = maximum(G.L) + 1
            J         = zeros(Int, mL, $D)
            B         = ones(mL, $D)
            dM = map(xij -> M(level(xij)), G.grid)
            yi = zeros(nA)

            @threadsfixed [J, B, yi] for i = 1:nx
                $(b[2])
                coverloop

                i1, i2 = UInt32(1), UInt32(1)
                iswitch = false

                @inbounds for ii = nc2
                    if G.adapt.overlap[i1][ii] && G.adapt.overlap[i2][ii]
                        b = 1.0
                        for d = 1:$D
                            b *= $(b[3])(x[i, d], G.grid[ii, d], dM[ii, d])
                            b == 0 && break
                        end
                        if b > 0
                            (yi += w[ii, :] * b)
                            iswitch = !iswitch
                            iswitch ? (i1 = ii::UInt32) : (i2 = ii::UInt32)
                        end
                    end
                end
                y[i, :] = yi
            end
            y
        end)

        subs!(f, :coverloop => coverloop)
        eval(current_module(), f)
    end
end
