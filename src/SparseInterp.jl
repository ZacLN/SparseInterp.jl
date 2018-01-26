__precompile__()
module SparseInterp
using Base.Threads
include("utils.jl")
include("univariategrids.jl")
include("grids.jl")
include("jlfuncs.jl")
include("adapt.jl")
include("api.jl")

export Linear, Quadratic,
       getW,
       interpolate,
       level,
       grow!,
       @threads,
       NGrid
end
