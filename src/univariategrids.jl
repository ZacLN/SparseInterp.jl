@inline M(i::T) where {T}           = (i == 1) ? T(1) : T(2)^(i - T(1)) + T(1)
@inline cc_iM(i::Int)           = (i == 1) ? 1 : Int(log2(i - 1) + 1)
@inline dM(i::T) where {T}           = (i < T(3)) ? T(2^(i - 1)) : T(2^(i - 2))
@inline cc_g(i::Int, j::Int)     = (i == 1) ?  0.5  : (j - 1) / (M(i) - 1.0)
@inline cc_g(i::Int)            = Float64[cc_g(i, j) for j = 1:M(i)]
@inline cc_dg(i::Int, j::Int)    = i == 2 ? [0.0, 1.0][j] : cc_g(i, 2j)
@inline cc_dg(i::Int)           = Float64[cc_dg(i, j) for j = 1:dM(i)]

@inline function cc_bf_l(x::Float64, xij::Float64, mi)
    if (mi == 1)
        return 1.0
    end
    dx = (1.0 - (mi - 1.0) * abs(x - xij))
    return (dx > 0.0) ? dx : 0.0
end

cc_bf_l(x::Float64, xij::Float64) = cc_bf_l(x, xij, M(level(xij)))
cc_bf_l(x::Vector{Float64}, xij::Vector{Float64}) = prod(map(cc_bf_l, x, xij))

@inline function cc_bf_q(x::Float64, xij::Float64, mi)
	if (mi == 1)
		return 1.0
	end
	dx = 1.0 - ((mi - 1.0) * (x - xij))^2
	return (dx > 0.0) ?  dx : 0.0
end

cc_bf_q(x::Float64, xij::Float64) = cc_bf_q(x, xij, cc_M(level(xij)))
cc_bf_q(x::Vector{Float64}, xij::Vector{Float64}) = map(cc_bf_q, x, xij)

function cc_simpsonsw(i::Int, j::Int)
    @assert 1 ≤ j ≤ M(i)
    if i == 1
        return 1.0
    end
    w = 1 / (M(i) - 1) / 3
    if j == 1 || j == M(i)
        return w
    else
        w *= mod(j, 2) == 0 ? 4 : 2
    end
    return w
end

cc_simpsonsw(i::Int) = Float64[cc_simpsonsw(i, j) for j = 1:M(i)]

function cc_dsimpsonsw(i::Int)
    if i == 1
        return cc_simpsonsw(i)
    elseif i == 2
        return cc_simpsonsw(i) - [0;1.0;0]
    else
        return cc_simpsonsw(i) - vec(Float64[cc_simpsonsw(i - 1)'zeros(1, M(i - 1))])[1:M(i - 1) - 1]
    end
end
