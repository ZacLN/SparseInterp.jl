function nXtoU(X::Array{Float64, 2}, bnds::Array{Float64, 2})
	N, R = size(X)
	out = Array{Float64}(size(X))
	for r = 1:R
		lb = bnds[1, r]
		ub = bnds[2, r]
		@inbounds for i = 1:N
			out[i, r] = (X[i, r] - lb) / (ub - lb)
		end
	end
	return out
end

function nUtoX(X::Array{Float64, 2}, bnds::Array{Float64, 2})
	N, R = size(X)
	out = Array{Float64}(size(X))
	for r = 1:R
		lb = bnds[1, r]
		ub = bnds[2, r]
		@inbounds for i = 1:N
			out[i, r] = X[i, r] * (ub - lb) + lb
		end
	end
	return out
end



"""
    subs!(x,list)

Substitute variables in an expression.
"""
function subs!(x::Expr, list::Dict)
    for i = 1:length(x.args)
        if in(x.args[i], keys(list))
            x.args[i] = list[x.args[i]]
        elseif isa(x.args[i], Expr)
            subs!(x.args[i], list)
        end
    end
end

subs!(x::Expr, p::Pair) = subs!(x, Dict(p))

function remlineinfo!(x)
    if isa(x, Expr)
        id = find(map(y -> isa(y, Expr) && y.head == :line, x.args))
        deleteat!(x.args, id)
        for j in x.args
            remlineinfo!(j)
        end
    end
end

# This alters Base.Threads.@threads to allow for individual copies of
# variables across a parallel loop:
# @threads [a,b...] for i = 1:...
# creates copies of a,b... before the parallel loop starts thus avoiding
# the need for multiple allocations
function _threadsforfixed(fixed, iter, lbody)
    fun = gensym("_threadsfor")
    lidx = iter.args[1]         # index
    range = iter.args[2]

    svars = map(gensym, fixed.args)
    newvar = quote end
    for i = 1:length(fixed.args)
        push!(newvar.args, :($(svars[i]) = copy($(fixed.args[i]))))
    end
    subs!(lbody, Dict(zip(fixed.args, svars)))

    quote
        function $fun()
            tid = threadid()
            r = $(esc(range))
            # divide loop iterations among threads
            len, rem = divrem(length(r), nthreads())
            # not enough iterations for all the threads?
            if len == 0
                if tid > rem
                    return
                end
                len, rem = 1, 0
            end
            # compute this thread's iterations
            f = 1 + ((tid - 1) * len)
            l = f + len - 1
            # distribute remaining iterations evenly
            if rem > 0
                if tid <= rem
                    f = f + (tid - 1)
                    l = l + tid
                else
                    f = f + rem
                    l = l + rem
                end
            end
            # run this thread's iterations

            $(esc(quote
            $newvar
        end))
            for i = f:l
                local $(esc(lidx)) = Base.unsafe_getindex(r, i)
                $(esc(lbody))
            end
        end
        ccall(:jl_threading_run, Void, (Any,), Core.svec($fun))
    end
end


# import Base.Threads.@threads
macro threadsfixed(args...)
    na = length(args)
    if na == 1
        ex = args[1]
        if !isa(ex, Expr)
            throw(ArgumentError("need an expression argument to @threads"))
        end
        if ex.head === :for
            return Base.Threads._threadsfor(ex.args[1], ex.args[2])
        else
            throw(ArgumentError("unrecognized argument to @threads"))
        end
    elseif na == 2
        fixed = args[1]
        ex = args[2]
        if ex.head === :for
            return _threadsforfixed(fixed, ex.args[1], ex.args[2])
        else
            throw(ArgumentError("unrecognized argument to @threads"))
        end
    else
        throw(ArgumentError("wrong number of arguments (1 or 2)"))
    end
end



"""
    comb(D,Q)

Generates a vector of  vectors containing D integers
that sum up to Q. The result is returned in decreasing
lexicographic order.
"""
function comb(D::Int, Q::Int)
	D == Q ? (return Vector{Int}[ones(Int, D)]) : nothing
	L = Q - D + 1
	out = Array{Vector{Int}}(binomial(Q - 1, D - 1))
	tL1 = ones(Int, D)
	tL2 = L * ones(Int, D)
	p = 1
	cnt = 0

	while tL1[D] < L
		tL1[p] += 1
		if tL1[p] > tL2[p]
			tL1[p] = 1
			p += 1
		else
			for i = 1:p - 1
				tL2[i] = tL2[p] - tL1[p] + 1end
			p = 1
			tL1[1] = tL2[1]
			cnt += 1
			out[cnt] = copy(tL1)
		end
	end
	return out
end

"""
    ndgrid(v...)

Computes the tensor grid of d input vectors and returns d
d-dimensional arrays.
"""
ndgrid(v::AbstractVector) = (copy(v),)

function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where {T}
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repmat(v1, 1, n), repmat(v2, m, 1))
end

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j - 1, snext), s) + 1]
    end
end

function ndgrid(vs::AbstractVector{T}...) where {T}
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i -> Array{T}(sz), n)
    s = 1
    for i = 1:n
        a = out[i]::Array
        v = vs[i]
        snext = s * size(a, i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
end

# Returns a tuple (l,j,dj) for any given point where
# l is the minimum level of grid to which it belongs,
# j is the index in this grid and dj the index in the
# difference grid.
function position(x::Float64)
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
    j = 0
    dj = 0
    if l == 1
        j = 1
        dj = 1
    elseif l == 2
        if x == 0.0
            j = 1
            dj = 1
        elseif x == 1.0
            j = 3
            dj = 2
        end
    else
        j = Int(div(x, 1 / 2^(l - 1)) + 1)
        dj = div(j - 1, 2) + 1
    end
    return l, j, dj
end

findrow(X::Array{T, 2}, x::Vector{T}) where {T} = findfirst(mapslices(hash, X, 2), hash(x))
