module EmbossPDE

import DoubleFloats
import GenericLinearAlgebra
import IntervalSets
import LinearAlgebra
import Polynomials
import Requires

using DoubleFloats: Double64
using IntervalSets: endpoints, (..)
using LinearAlgebra: dot, svd!, qr!, diag
using Polynomials: Polynomial, ChebyshevT

import Base: +, -, *

export operators, variables, solve, (..)

struct PDESolution{FloatType, Iter, CacheNamedTuple} <: Function
    n::Int
    c::Vector{FloatType}
    dom::Iter
    cache::CacheNamedTuple
end
coeffs(u::PDESolution) = u.c
domain(u::PDESolution) = domget.(u.dom)
domget(d::Base.Fix2) = (d.f, d.x)   # convert ≤(blah) to (≤, blah)
domget(d) = d

struct SubstitutionOperator{FloatType}
    n::Int
    nodes::Vector{FloatType}
end
(B::SubstitutionOperator)(s) = subs(s, B.n, B.nodes; prefun=identity) # y => f(x)
(B::SubstitutionOperator)(s::Tuple) = begin # y => f(x), a..b
    a,b = endpoints(last(s))
    prefun = node -> (node+1)/2 * (b-a) + a  # map (-1,1) to (a,b)
    subs(first(s), B.n, B.nodes; prefun)
end

struct MatrixAndFunction{FloatType, FunctionType}
    A::Matrix{FloatType}
    f::FunctionType
    var::Symbol  # the variable that f is a function of
end
function (Base.:*)(mf::MatrixAndFunction, other::Union{Number, AbstractMatrix}) # e.g. B(blah)*∂x
    MatrixAndFunction(mf.A*other, mf.f, mf.var)
end
function (Base.:-)(mf::MatrixAndFunction, mf2::MatrixAndFunction) # e.g. B(blah) - B(bleh)
    @assert mf.f == mf2.f
    @assert mf.var == mf2.var
    MatrixAndFunction(mf.A-mf2.A, mf.f, mf.var)
end
function (Base.:-)(mf::MatrixAndFunction) # e.g. -B(blah)
    MatrixAndFunction(-mf.A, mf.f, mf.var)
end
function (Base.:+)(mf::MatrixAndFunction, mf2::MatrixAndFunction) # e.g. B(blah) + B(bleh)  (but why?)
    @assert mf.f == mf2.f
    @assert mf.var == mf2.var
    MatrixAndFunction(mf.A+mf2.A, mf.f, mf.var)
end

chebT(FloatType, n) = ChebyshevT(FloatType.(LinearAlgebra.I(n+1)[n+1,:]))
chebT(FloatType, n, symbol) = ChebyshevT(FloatType.(LinearAlgebra.I(n+1)[n+1,:]), symbol)

function diffx(FloatType, n; degree::Integer=1)
    D1d = D1dmat(FloatType, n)^degree
    D = kron(D1d,LinearAlgebra.I(n+1))
    idx = (i+j for i=0:n for j=0:n) .≤ n
    D[idx, idx]
end

function diffy(FloatType, n; degree::Integer=1)
    D1d = D1dmat(FloatType, n)^degree
    D = kron(LinearAlgebra.I(n+1),D1d)
    idx = (i+j for i=0:n for j=0:n) .≤ n
    D[idx, idx]
end

function mulx(FloatType, n)
    X1d = X1dmat(FloatType, n)
    X = kron(X1d,LinearAlgebra.I(n+1))
    idx = (i+j for i=0:n for j=0:n) .≤ n
    X[idx, idx]
end

function muly(FloatType, n)
    Y1d = X1dmat(FloatType, n)
    Y = kron(LinearAlgebra.I(n+1),Y1d)
    idx = (i+j for i=0:n for j=0:n) .≤ n
    Y[idx, idx]
end

function mulfx(FloatType, f, n, nodes)
    fX1d = fX1dmat(FloatType, f, n, nodes)
    fX = kron(fX1d,LinearAlgebra.I(n+1))
    idx = (i+j for i=0:n for j=0:n) .≤ n
    fX[idx, idx]
end

function mulfy(FloatType, f, n, nodes)
    fY1d = fX1dmat(FloatType, f, n, nodes)
    fY = kron(LinearAlgebra.I(n+1),fY1d)
    idx = (i+j for i=0:n for j=0:n) .≤ n
    fY[idx, idx]
end

function D1dmat(FloatType, n)
    D = zeros(FloatType, n+1,n+1)
    for j = 1:n+1
        p = chebT(FloatType, j-1)
        cd2p = Polynomials.coeffs(Polynomials.derivative(p))
        D[1:length(cd2p),j] = cd2p
    end
    round.(D)   # exact values are integers
end

function X1dmat(FloatType, n)
    X = zeros(FloatType, n+1,n+1)
    for j = 1:n+1
        p = chebT(FloatType, j-1)
        cxp = Polynomials.coeffs(ChebyshevT((0,1))*p)
        len = min(n+1, length(cxp))
        X[1:len,j] = cxp[1:len]
    end
    X
end

function fX1dmat(FloatType, f, n, nodes)
    T(n,x) = cos(n*acos(x))
    W = reduce(hcat, T.(i,nodes) for i=0:n)
    fX1d = W'*(W.*f.(nodes))/n
    fX1d[1,:] ./= 2
    fX1d
end

function project(FloatType, f::Function, n, var, prefun)
    np = 2*n
    prenodes = [(2 * FloatType(k) - 1) * FloatType(π) / (2 * FloatType(np)) for k = np:-1:1]
    nodes = cos.(prenodes)
    project(f, n, nodes, var, prefun)
end

# f could be a function of x or y or x and y
function project(f::Function, n::Integer, nodes, var, prefun)
    @debug "project(function)"
    T(n,x) = cos(n*acos(x))
    W = reduce(hcat, T.(i,nodes) for i=0:n)
    rows = 1 .+ [0;cumsum(n+1:-1:2)]
    if var == :x₁₂ # f(x,y)
        F = [f(xᵢ, yⱼ) for xᵢ in nodes, yⱼ in nodes]
        c = [4*sum((W[:,i+1] * W[:,j+1]')/length(nodes)^2 .* F) for i=0:n for j=0:n if i+j ≤ n] # double integral over x,y
        c[1:n+1] /= 2   # for y
        c[rows] /= 2    # for x
    elseif var == :x₁ # f(x), cheap way
        cx = 2*W'*f.(prefun.(nodes))/length(nodes)
        cx[1] /= 2
        c = zeros(eltype(cx), (n+2)*(n+1)÷2)
        c[rows] .= cx
    elseif var == :x₂ # f(y), cheap way
        cy = 2*W'*f.(prefun.(nodes))/length(nodes)
        cy[1] /= 2
        c = zeros(eltype(cy), (n+2)*(n+1)÷2)
        c[1:n+1] .= cy
    else
        error("Internal error: confused about variable of substitution.")
    end
    c
end

function project(FloatType, p::Polynomial{T, :x₁}, n, var, prefun) where T
    @debug "project(poly(x)) fast path"
    if var == :x₂  # but we can see that the polynomial is in terms of x₁
        error("The function should be of :x₂")
    end
    @assert length(p) ≤ n
    rows = 1 .+ [0;cumsum(n+1:-1:2)]
    N = (n+1)*(n+2)÷2
    pp = p(prefun(Polynomial(:x₁)))
    pc = convert(ChebyshevT{FloatType, :x₁}, pp)
    c = Polynomials.coeffs(pc)
    C = zeros(eltype(c), N)
    len = min(length(c),n+1)
    C[rows[1:len]] .= c[1:len]
    C
end

function project(FloatType, p::Polynomial{T, :x₂}, n, var, prefun) where T
    @debug "project(poly(y)) fast path"
    if var == :x₁  # but we can see that the polynomial is in terms of x₂
        error("The function should be of :x₁")
    end
    @assert length(p) ≤ n
    N = (n+1)*(n+2)÷2
    pp = p(prefun(Polynomial(:x₂)))
    pc = convert(ChebyshevT{FloatType, :x₂}, pp)
    c = Polynomials.coeffs(pc)
    C = zeros(eltype(c), N)
    len = min(length(c),n+1)
    C[1:len] .= c[1:len]
    C
end

function project(FloatType, val::Number, n, var, prefun)
    @debug "project(number) fast path"
    N = (n+1)*(n+2)÷2
    C = zeros(FloatType, N)
    C[1] = val
    C
end

function subs(s::Pair{Polynomial{T, :x₁}, F}, n, nodes; prefun) where {T,F}
    if s[1] == Polynomials.variable(:x₁)  # and not, say, x^2 or something
        subs(:x₁=>s[2], n, nodes; prefun)
    else
        error("You can only substitute for x₁ or x₂.")
    end
end

function subs(s::Pair{Polynomial{T, :x₂}, F}, n, nodes; prefun) where {T,F}
    if s[1] == Polynomials.variable(:x₂)  # and not, say, y^2 or something
        subs(:x₂=>s[2], n, nodes; prefun)
    else
        error("You can only substitute for x₁ or x₂.")
    end
end

check_var_constency(var::Symbol, ::Polynomial{T, :x₁}) where T = var == :x₂
check_var_constency(var::Symbol, ::Polynomial{T, :x₂}) where T = var == :x₁
check_var_constency(var::Symbol, _) = true  # for all we know

function subs(s::Pair{Symbol, F}, n, nodes; prefun) where F
    if !check_var_constency(s[1], s[2])
        error("Invalid variable substitution: same variable on either side.")
    end
    if s[1] == :x₁
        subsx(s[2], n, nodes; prefun)
    elseif s[1] == :x₂
        subsy(s[2], n, nodes; prefun)
    else
        error("Symbol must be :x₁ or :x₂.")
    end
end

function subsx(funy, n, nodes; prefun)
    @debug "subsx(funy)"
    N = (n+1)*(n+2)÷2
    FloatType = typeof(nodes[1]*prefun(nodes[1])*funy(prefun(nodes[1])))
    A = zeros(FloatType, N,N)
    T(n,x) = cos(n*acos(x))
    W = reduce(hcat, T.(i,nodes) for i=0:n)
    V = reduce(hcat, T.(i,prefun.(nodes)) for i=0:n)
    Y = reduce(hcat, T.(i,funy.(prefun.(nodes))) for i=0:n)
    # S = reduce(hcat, Y[:,i+1] .* V[:,j+1] for i=0:n for j=0:n if i+j ≤ n)  # too slow
    idx = [(i,j) for i=0:n for j=0:n if i+j ≤ n]
    S = zeros(FloatType, length(nodes), length(idx))
    Threads.@threads for k in eachindex(idx)
        (i,j) = idx[k]
        @views S[:,k] .= Y[:,i+1] .* V[:,j+1]
    end
    A[1:n+1, :] .= 2*W'*S/length(nodes)
    A[1,:] ./= 2
    MatrixAndFunction(A, prefun, :x₂)
end

function subsy(funx, n, nodes; prefun)
    @debug "subsy(funx)"
    N = (n+1)*(n+2)÷2
    FloatType = typeof(nodes[1]*prefun(nodes[1])*funx(prefun(nodes[1])))
    A = zeros(FloatType, N,N)
    rows = 1 .+ [0;cumsum(n+1:-1:2)]
    T(n,x) = cos(n*acos(x))
    W = reduce(hcat, T.(i,nodes) for i=0:n)
    V = reduce(hcat, T.(i,prefun.(nodes)) for i=0:n)
    Y = reduce(hcat, T.(i,funx.(prefun.(nodes))) for i=0:n)
    # S = reduce(hcat, V[:,i+1] .* Y[:,j+1] for i=0:n for j=0:n if i+j ≤ n) # too slow
    idx = [(i,j) for i=0:n for j=0:n if i+j ≤ n]
    S = zeros(FloatType, length(nodes), length(idx))
    Threads.@threads for k in eachindex(idx)
        (i,j) = idx[k]
        @views S[:,k] .= V[:,i+1] .* Y[:,j+1]
    end
    A[rows, :] .= 2*W'*S/length(nodes)
    A[1,:] ./= 2
    MatrixAndFunction(A, prefun, :x₁)
end

function subsx(val::Number, n, nodes; prefun)
    @debug "subsx(number) fast path"
    FloatType = eltype(nodes)
    floatval = FloatType(val)
    subsx(y->floatval, n, nodes; prefun)
end

function subsy(val::Number, n, nodes; prefun)
    @debug "subsy(number) fast path"
    FloatType = eltype(nodes)
    floatval = FloatType(val)
    subsy(x->floatval, n, nodes; prefun)
end

function operators(FloatType, n; maxdegree::Integer=1)
    np = 2*n
    prenodes = [(2 * FloatType(k) - 1) * FloatType(π) / (2 * FloatType(np)) for k = np:-1:1]
    nodes = cos.(prenodes)
    Dx = diffx(FloatType, n)
    Dy = diffy(FloatType, n)

    if maxdegree > 1
        Dx = tuple(Dx, (diffx(FloatType, n; degree=k) for k=2:maxdegree)...)
        Dy = tuple(Dy, (diffy(FloatType, n; degree=k) for k=2:maxdegree)...)
    end

    B = SubstitutionOperator{FloatType}(n, nodes)

    X = mulx(FloatType, n)
    Y = muly(FloatType, n)

    Fx = f -> mulfx(FloatType, f, n, nodes)
    Fy = f -> mulfy(FloatType, f, n, nodes)

    Dx, Dy, B, X, Y, Fx, Fy
end

operators(n; maxdegree=1) = operators(Float64, n; maxdegree)

variables() = Polynomials.variable.((:x₁, :x₂))

# TODO: think more about how best to do this type promotion stuff
getrhstype(T, x::Number, var, prefun) = typeof(x)
getrhstype(T, p::Polynomial, var, prefun) = T
function getrhstype(T, f::Function, var, prefun)
    if var == :x₁₂
        typeof(f(zero(T),zero(T)))
    else
        typeof(f(prefun(zero(T))))
    end
end

function isemptyrow(row, b)
    iszero(b) && all(iszero, row)
end

function isemptyrow(row, b, tol)
    abs(b) ≤ tol && maximum(abs, row) ≤ tol
end

# a bit hacky
stripeqns(A) = A, :x₁₂, nothing
stripeqns(mfp::Pair{MatrixAndFunction{M,F}, T}) where {M,F,T} = first(mfp).A => last(mfp), first(mfp).var, first(mfp).f

function assemble(raweqns...; remove_empty_rows=true)

    # the code just sort of evolved to this point, which looks a bit hacky
    eqns_and_prefuns = stripeqns.(raweqns)
    eqns = getindex.(eqns_and_prefuns, 1)
    vars = getindex.(eqns_and_prefuns, 2)
    prefuns = getindex.(eqns_and_prefuns, 3)

    N = size(first(eqns[1]), 1)  # all matrices are the same size
    n = round(Int, sqrt(1/4 + 2*N) - 3/2)
    neq = length(eqns)

    # TODO: think more about how best to do this type promotion stuff
    MatrixFloatType = typeof(sum(zero.(eltype.(first.(eqns)))))
    RHSTypes = getrhstype.(MatrixFloatType, last.(eqns), vars, prefuns)
    RHSFloatType = typeof(zero(MatrixFloatType) + sum(zero.(RHSTypes)))

    # Assemble right hand side first, because...
    B = zeros(RHSFloatType, N, length(eqns))
    for (j,fun) in pairs(last.(eqns))
        B[:,j] = project(RHSFloatType, fun, n, vars[j], prefuns[j])
    end

    # ...we might want to omit rows from [A b] that are entirely zeros
    keeprow = trues(N, neq)
    Ascale = maximum(M->maximum(abs, M), first.(eqns))
    bscale = maximum(abs, B)
    problemscale = max(Ascale, bscale)
    @debug problemscale
    if remove_empty_rows
        cutoff = eps(MatrixFloatType) * problemscale  # is this a good cutoff?
        for j = 1:neq
            Aj = first(eqns[j])
            Threads.@threads for i = 1:N
                @views if isemptyrow(Aj[i,:], B[i,j], cutoff)
                    keeprow[i,j] = false
                end
            end
        end
    end

    # Now copy across the rows that we want to keep.
    nkeepeq = sum(keeprow)
    A = Matrix{MatrixFloatType}(undef, nkeepeq, N)
    b = Vector{RHSFloatType}(undef, nkeepeq)
    for j = 1:neq
        Aj = first(eqns[j])
        blockidxs = [0; cumsum(sum(keeprow, dims=1)[:])]
        @assert blockidxs[end] == nkeepeq
        @views A[blockidxs[j]+1:blockidxs[j+1], :] = Aj[keeprow[:,j], :]
        b[blockidxs[j]+1:blockidxs[j+1]] = B[keeprow[:,j], j]
    end

    # Check we didn't stuff up the indexing
    # @assert isequal(reduce(vcat, M for M in first.(eqns))[keeprow[:],:], A)
    # @assert isequal(B[keeprow[:]], b)

A, b, keeprow[:]

end

function factorise!(A; method)
    if method == :qr
        F = qr!(A, LinearAlgebra.ColumnNorm())
        svals = abs.(diag(F.R))  # sort of not really
    elseif method == :svd
        F = svd!(A)
        svals = F.S
    else
        error("Bad method.")
    end
    F, svals
end

function linsolve(F, b; method, rtol)
    if method == :qr
        qrsolve(F, b; rtol)
    elseif method == :svd
        svdsolve(F, b; rtol)
    else
        error("Bad method.")
    end
end

function qrsolve(F, b; rtol)
    # truncated pivoted QR solve; Julia Base doesn't support custom rtol
    m,n = size(F)
    @assert m ≥ n
    k = searchsortedlast(abs.(diag(F.R)), abs(F.R[1])*rtol, rev=true)
    y = (F.Q' * b)[1:k]
    z = @views [F.R[1:k,1:k] \ y; zeros(n-k)]
    z[invperm(F.p)]
end

function svdsolve(F, b; rtol)
    # truncated SVD solve; Julia Base doesn't support custom rtol
    k = searchsortedlast(F.S, F.S[1]*rtol, rev=true)
    @views F.Vt[1:k,:]' * (F.S[1:k] .\ (F.U[:,1:k]' * b))
end

extended_precision(x::Float32) = Float64(x)
extended_precision(x::Float64) = Double64(x)
extended_precision(x::Double64) = BigFloat(x)
extended_precision(x) = one(extended_precision(x.value)) * x  # for ForwardDiff

function solve(raweqns...; domain, method=:qr, cutoff=eps, remove_empty_rows=true, refine=false)

    if method ∉ (:qr, :svd)
        error("Invalid method.  Must be :qr or :svd")
    end

    A,b,idx = assemble(raweqns...; remove_empty_rows)

    eqns = first.(stripeqns.(raweqns))

    FloatType = eltype(A)
    rtol = cutoff == eps ? eps(FloatType) : cutoff

    F, svals = factorise!(A; method)
    c = linsolve(F, b; method, rtol)

    CoefficientType = eltype(c)

    if refine # iterative refinement in extended precision
        # we overwrote A to save memory, but we can still get the residual from the original inputs
        bfit_extended_precision = reduce(vcat, (M*extended_precision.(c) for M in first.(eqns)))
        r = extended_precision.(b) - bfit_extended_precision[idx]
        δc = linsolve(F, r; method, rtol) # n.b. rtol is _not_ the extended precision rtol
        c = CoefficientType.(c + δc)
    end

    # compute fitted b based on the full inputs as a further check we didn't stuff up
    bfit = reduce(vcat, (M*c for M in first.(eqns))) # no omitting zero rows here
    # @assert all(iszero, bfit[.!idx])   # TODO: re-enable this check with a suitable tolerance

    N = sum(size(M,1) for M in first.(eqns))
    maxresidual = maximum(abs, b - bfit[idx]) / maximum(abs, b) # relative to scale of b
    @info "Solve stats" N size(A) cond(A)=svals[1]/svals[end] maxresidual
    cache = (; svals, maxresidual, cutoff=svals[1]*rtol)

    reconstitute(c, domain; cache)
end

function reconstitute(c, domain; cache)
    N = length(c)
    n = round(Int, sqrt(1/4 + 2*N) - 3/2)
    PDESolution(n, c, domain, cache)
end

function check_condition(op, var1, var2, fun, range)
    var2 ∉ range || op(var1, fun(var2))
end

function check_condition(op, var1, var2, value::Number, range)
    var2 ∉ range || op(var1, value)
end

function raweval(u, x, y)
    n = u.n
    dot(u.c, (cos(i*acos(x))*cos(j*acos(y)) for i=0:n for j=0:n if i+j ≤ n))
end

function destructure_boundary(boundary)
    op, varfunrange = boundary
    if varfunrange isa Tuple
        var, fun = first(varfunrange)
        range = last(varfunrange)
    else
        var, fun = varfunrange
        range = -1..1
    end
    op, var, fun, range
end

function nanwrapeval(f, x, y, domain)
    for boundary in domain
        op, var, fun, range = destructure_boundary(boundary)
        if var == Polynomials.variable(:x₁)
            if !check_condition(op, x, y, fun, range)
                return NaN
            end
        elseif var == Polynomials.variable(:x₂)
            if !check_condition(op, y, x, fun, range)
                return NaN
            end
        else
            error("Boundary must be in terms of a single variable only.")
        end
    end
    f(x,y)
end

function (u::PDESolution{FloatType})(x,y; mask=true)::FloatType where FloatType
    if mask
        nanwrapeval((x,y)->raweval(u, x, y), x, y, domain(u))
    else
        raweval(u, x, y)
    end
end

function __init__()
    Requires.@require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" Makie.plot(
        u::PDESolution; divisions=256, levels=20, size=(1200,600),
        aspect=1, axis=(;aspect), diagnostics=true, mask=true
    ) = 
    begin
        gs = divisions
        FloatType = eltype(u.c)
        grid = range(-1,1,length=gs)
        Z = zeros(gs,gs)
        Threads.@threads for i=1:gs for j=1:gs
            Z[i,j] = u(grid[i], grid[j]; mask)
        end end
        fig = Makie.Figure(;size)
        ax, hm = Makie.heatmap(fig[1:2,1:2], grid, grid, Z; axis)
        ax.title = "max |rₖ| / max |bₖ| = " * string(Float32(u.cache.maxresidual))
        Makie.contour!(fig[1:2,1:2], grid,grid,Z; color=:white, levels)
        Makie.Colorbar(fig[1:2,3], hm)
        if diagnostics
            N = length(u.c)
            n = round(Int, sqrt(1/4 + 2*N) - 3/2)
            ax3, hm3 = Makie.scatter(fig[1,4], abs.(u.c).+eps(FloatType), markersize=6, axis=(xlabel="𝑘", ylabel="|𝐶ₖ|", yscale=log10))
            ax3.title = "n = " * string(n) * ";  " * string(N) * " terms"
            ax2, hm2 = Makie.scatter(fig[2,4], u.cache.svals.+eps(FloatType), markersize=6, axis=(xlabel="𝑘", ylabel="𝑆ₖ", yscale=log10))
            Makie.lines!(fig[2,4], [0,N], u.cache.cutoff*[1,1]; color=:red)
            ax2.title = "cond(A) ≈ " * string(Float32(u.cache.svals[1] / u.cache.svals[end]))
        end
        fig
    end
end

end # module EmbossPDE
