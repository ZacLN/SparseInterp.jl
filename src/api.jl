(G::NGrid)(A::Vector{Float64}, x::Array{Float64, 2}, y = zeros(size(x, 1))) = jl_ainterp(G, A, x, y)

(G::NGrid)(A::Array{Float64, 2}, x::Array{Float64, 2}, y = zeros(size(x, 1), size(A, 2))) = jl_ainterp(G, A, x, y)
