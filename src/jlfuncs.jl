"""
    buildW(G,hind,hcovers,pr)

Precomputes
"""
function buildW(G::NGrid{D, BF}, hind, hcovers, pr = 1:length(G)) where {D,BF}
    bf = (BF == Linear ? cc_bf_l : cc_bf_q)
    nc 		= size(G.covers, 1)

    coversid = zeros(Int32, length(G))
    for i = 1:length(hind)
        coversid[i] = findfirst(hcovers, hind[i])
    end
    mL      = maximum(G.L)
    J 		= ones(Int, mL + 1, D)
    B 		= zeros(mL + 1, D)
    for i = pr
        G.IDs[i] = Int[]
        G.Bs[i] = Float64[]
        for d = 1:D
            for l = 1:mL + 1
                j 	= clamp(round(Int, G.grid[i, d] * (dM(l)) + 1 / 2), 1, dM(l))
                B[l, d] 	= bf(G.grid[i, d], cc_dg(l, j), Int16(M(l)))
                J[l, d]  = j
            end
        end
        for ii = 1:coversid[i]
            b  = B[G.covers[ii, D], D] * B[G.covers[ii, 1], 1]
            id1 = J[G.covers[ii, D], D] - 1
            for d = D - 1:-1:2
                b *= B[G.covers[ii, d], d]
                id1 = id1 * G.covers_dM[ii, d] + (J[G.covers[ii, d], d] - 1)
            end
            id1 = (J[G.covers[ii, 1], 1] - 1) + G.covers_dM[ii, 1] * id1 + 1 + G.covers_loc[ii] - 1
            if b > 0
                push!(G.IDs[i], id1)
                push!(G.Bs[i], id1 == i ? b : -b)
            end
        end
    end
    return
end


# This section constructs specialised functions for grids of dimensions 2-12
Lbj = quote
@inbounds for d in dr1
    J[2, d] = (x[i, d] > 0.5)
    B[2, d] = (1.0 - 2 * abs(x[i, d] - J[2, d]))
end
@inbounds for l = 3:mL
    dm = 2^(l - 2)
    m  = dm * 2
    for d in dr1
        j 	= clamp(round(Int, x[i, d] * dm + 1 / 2), 1, dm)
        xij = (2j - 1) / m
        dx = (1.0 - (m) * abs(x[i, d] - xij))
        B[l, d] 	= dx
        J[l, d]  = j - 1
    end
end
end

Qbj = quote
@inbounds for d in dr1
    J[2, d] = (x[i, d] > 0.5)
    B[2, d] = 1.0 - (2 * (x[i, d] - J[2, d]))^2
end
@inbounds for l = 3:mL
    dm = 2^(l - 2)
    m  = dm * 2
    for d in dr1
        j 	= clamp(round(Int, x[i, d] * dm + 1 / 2), 1, dm)
        xij = (2j - 1) / m
        dx = 1.0 - (m * (x[i, d] - xij))^2
        B[l, d]  = dx
        J[l, d]  = j - 1
    end
end
end



for b in [(Linear, Lbj), (Quadratic, Qbj)]
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

        f = :(function interpolate(G::NGrid{$D, $(b[1])}, A::Vector{Float64}, xi::Array{Float64, 2}, y = zeros(Float64, size(xi, 1)))
            w         = getW(G, A)
            x         = nXtoU(xi, G.bounds)
            nx        = size(x, 1)
            cr         = 1:size(G.covers, 1)
            dr1 = 1:$D
            mL      = maximum(G.L) + 1
            J         = zeros(Int, mL, $D)
            B         = ones(mL, $D)

            @threadsfixed [J, B] for i = 1:nx
                $(b[2])
                yi = 0.0
                @inbounds for ii in cr
                    $(Expr(:block,
                        :(b   = B[G.covers[ii, $D], $D] * B[G.covers[ii, 1], 1]),
                        (:(b  *= B[G.covers[ii, $d], $d]) for d in D - 1:-1:2)...,
                        :(id1 = J[G.covers[ii, $D], $D]),
                        (:(id1 = id1 * G.covers_dM[ii, $d] + J[G.covers[ii, $d], $d]) for d in D - 1:-1:2)...,
                        :(id1 = id1 * G.covers_dM[ii, 1] + G.covers_loc[ii] + J[G.covers[ii, 1], 1]),
                        :(yi += b * w[id1])
                        ))
                end
                y[i] = yi
            end
            y
        end)

        subs!(f, :coverloop => coverloop)
        eval(current_module(), f)
    end
end

mutable struct dimdef{n} end


for b in [(Linear, Lbj), (Quadratic, Qbj)]
    for D = 2:15
        for adim = 1:8
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
            for d = 1:adim
                push!(coverloop.args[4].args[2].args[2].args, :(yi[$d] += b * w[id1, $d]))
                push!(coverloop.args, :(@inbounds y[i, $d] = yi[$d]))
            end

            f = :(function interpolate(G::NGrid{$D, $(b[1])}, A::Array{Float64, 2}, xi::Array{Float64, 2}, nA::dimdef{$adim}, y::Array{Float64, 2})
                w         = getW(G, A)
                x         = nXtoU(xi, G.bounds)
                nx        = size(x, 1)

                nc         = 1:size(G.covers, 1)
                dr1 = 1:$D
                mL      = maximum(G.L) + 1
                J         = zeros(Int, mL, $D)
                B         = ones(mL, $D)
                yi        = zeros($adim)

                @threadsfixed [J, B, yi] for i = 1:nx
                    $(b[2])
                    coverloop
                end
                y
            end)

            subs!(f, :coverloop => coverloop)
            eval(current_module(), f)
        end
    end
end

interpolate(G::NGrid, A::Array{Float64, 2}, x::Array{Float64, 2}, y = zeros(size(x, 1), size(A, 2))) = interpolate(G, A, x, dimdef{size(A, 2)}(), y)


function getW(G::NGrid, A::Vector{Float64})
    w = copy(A)
    for i = 1:length(G)
        for l ∈ 1:length(G.Bs[i]) - 1
            @inbounds w[i] +=  G.Bs[i][l] * w[G.IDs[i][l]]
        end
    end
    return w
end

function getW(G::NGrid, A::Array{Float64, 2})
    nA = size(A, 2)
    w = copy(A)
    for i = 1:length(G)
        for l ∈ 1:length(G.Bs[i]) - 1
            for a = 1:nA
                @inbounds w[i, a] +=  G.Bs[i][l] * w[G.IDs[i][l], a]
            end
        end
    end
    return w
end
