using SparseInterp, Base.Test

G = NGrid([8,8,8,8,8])
A = rand(length(G))

@test all(isapprox.(G(A, G.grid), A))