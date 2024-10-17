module EmbossPDE

import DoubleFloats
import GenericLinearAlgebra
import IntervalSets
import LinearAlgebra
import Polynomials
import Requires

using DoubleFloats: Double64
using IntervalSets: endpoints, width, (..)
using LinearAlgebra: dot, svd!, qr!, diag, norm
using Polynomials: Polynomial, ChebyshevT

import Base: +, -, *

export operators, variables, solve, (..)

struct PDESolution{FloatType, Iter, CacheNamedTuple} <: Function
    n::Int
    c::Vector{FloatType}
    boundingbox::Tuple{IntervalSets.ClosedInterval{FloatType}, IntervalSets.ClosedInterval{FloatType}}
    dom::Iter
    cache::CacheNamedTuple
end
coeffs(u::PDESolution) = u.c
domain(u::PDESolution) = domget.(u.dom)
boundingbox(u::PDESolution) = u.boundingbox
domget(d::Base.Fix2) = (d.f, d.x)   # convert ≤(blah) to (≤, blah)
domget(d) = d

struct SubstitutionOperator{FloatType}
    n::Int
    nodes::Vector{FloatType}
    boundingbox::Tuple{IntervalSets.ClosedInterval{FloatType}, IntervalSets.ClosedInterval{FloatType}}
end
(B::SubstitutionOperator)(s) = subs(s, B.n, B.nodes; prefun=identity, boundingbox=B.boundingbox) # y => f(x)
function (B::SubstitutionOperator)(s::Tuple) # x => f(y), a..b  or  y => f(x), a..b
    var = first(first(s))
    if var == Polynomial(:x₁)
        aa, bb = endpoints(last(B.boundingbox))  # careful here, if x is the variable, then the bounds are for y
    elseif var == Polynomial(:x₂)
        aa, bb = endpoints(first(B.boundingbox))
    else
        error("Internal error: Confused by variable substitution.")
    end

    to_minus_one_one = var -> (bb + aa - 2*var)/(-bb + aa)

    a,b = to_minus_one_one.(endpoints(last(s)))
    prefun = node -> (node+1)/2 * (b-a) + a  # map (-1,1) to (a,b)
    subs(first(s), B.n, B.nodes; prefun, boundingbox=B.boundingbox)
end

struct MatrixAndFunction{FloatType, FunctionType}
    A::Matrix{FloatType}
    f::FunctionType
    var::Symbol  # the variable that f is a function of
    boundingbox::Tuple{IntervalSets.ClosedInterval{FloatType}, IntervalSets.ClosedInterval{FloatType}}
end
MatrixAndFunction(A, f, var, boundingbox) = MatrixAndFunction{eltype(A), typeof(f)}(A, f, var, boundingbox)

function (Base.:*)(mf::MatrixAndFunction, other::Union{Number, AbstractMatrix}) # e.g. B(blah)*∂x
    MatrixAndFunction(mf.A*other, mf.f, mf.var, mf.boundingbox)
end
function (Base.:-)(mf::MatrixAndFunction, mf2::MatrixAndFunction) # e.g. B(blah) - B(bleh)
    @assert mf.f == mf2.f
    @assert mf.var == mf2.var
    @assert mf.boundingbox == mf2.boundingbox
    MatrixAndFunction(mf.A-mf2.A, mf.f, mf.var, mf.boundingbox)
end
function (Base.:-)(mf::MatrixAndFunction) # e.g. -B(blah)
    MatrixAndFunction(-mf.A, mf.f, mf.var, mf.boundingbox)
end
function (Base.:+)(mf::MatrixAndFunction, mf2::MatrixAndFunction) # e.g. B(blah) + B(bleh)  (but why?)
    @assert mf.f == mf2.f
    @assert mf.var == mf2.var
    @assert mf.boundingbox == mf2.boundingbox
    MatrixAndFunction(mf.A+mf2.A, mf.f, mf.var, mf.boundingbox)
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

function mulfx(FloatType, f, n, nodes, boundingbox)
    fX1d = fX1dmat(FloatType, f, n, nodes, boundingbox)
    fX = kron(fX1d,LinearAlgebra.I(n+1))
    idx = (i+j for i=0:n for j=0:n) .≤ n
    fX[idx, idx]
end

function mulfy(FloatType, f, n, nodes, boundingbox)
    fY1d = fX1dmat(FloatType, f, n, nodes, boundingbox)
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

function fX1dmat(FloatType, f, n, nodes, boundingbox)
    ax, bx = endpoints(first(boundingbox))
    x_from_minus_one_one = x -> ax - (1 + x)*(-bx/2 + ax/2)
    x_to_minus_one_one = x -> (bx + ax - 2*x)/(-bx + ax)
    T(n,x) = cos(n*acos(x))
    W = reduce(hcat, T.(i,nodes) for i=0:n)
    fX1d = W'*(W.*f.(x_from_minus_one_one.(nodes)))/n
    fX1d[1,:] ./= 2
    fX1d
end

function project(FloatType, f::Function, n, var, prefun, boundingbox)
    np = 2*n
    prenodes = [(2 * FloatType(k) - 1) * FloatType(π) / (2 * FloatType(np)) for k = np:-1:1]
    nodes = cos.(prenodes)
    project(f, n, nodes, var, prefun, boundingbox)
end

# f could be a function of x or y or x and y
function project(f::Function, n::Integer, nodes, var, prefun, boundingbox)
    @debug "project(function)"
    T(n,x) = cos(n*acos(x))
    W = reduce(hcat, T.(i,nodes) for i=0:n)
    rows = 1 .+ [0;cumsum(n+1:-1:2)]

    ax, bx = endpoints(first(boundingbox))
    x_from_minus_one_one = x -> ax - (1 + x)*(-bx/2 + ax/2)
    x_to_minus_one_one = x -> (bx + ax - 2*x)/(-bx + ax)

    ay, by = endpoints(last(boundingbox))
    y_from_minus_one_one = y -> ay - (1 + y)*(-by/2 + ay/2)
    y_to_minus_one_one = y -> (by + ay - 2*y)/(-by + ay)

    if var == :x₁₂ # f(x,y)
        F = [f(x_from_minus_one_one(xᵢ), y_from_minus_one_one(yⱼ)) for xᵢ in nodes, yⱼ in nodes]
        c = [4*sum((W[:,i+1] * W[:,j+1]')/length(nodes)^2 .* F) for i=0:n for j=0:n if i+j ≤ n] # double integral over x,y
        c[1:n+1] /= 2   # for y
        c[rows] /= 2    # for x
    elseif var == :x₁ # f(x), cheap way
        cx = 2*W'*f.(x_from_minus_one_one.(prefun.(nodes)))/length(nodes)
        cx[1] /= 2
        c = zeros(eltype(cx), (n+2)*(n+1)÷2)
        c[rows] .= cx
    elseif var == :x₂ # f(y), cheap way
        cy = 2*W'*f.(y_from_minus_one_one.(prefun.(nodes)))/length(nodes)
        cy[1] /= 2
        c = zeros(eltype(cy), (n+2)*(n+1)÷2)
        c[1:n+1] .= cy
    else
        error("Internal error: confused about variable of substitution.")
    end
    c
end

function project(FloatType, p::Polynomial{T, :x₁}, n, var, prefun, boundingbox) where T
    @debug "project(poly(x)) fast path"
    if var == :x₂  # but we can see that the polynomial is in terms of x₁
        error("The function should be of :x₂")
    end
    @assert length(p) ≤ n

    ax, bx = endpoints(first(boundingbox))
    x_from_minus_one_one = x -> ax - (1 + x)*(-bx/2 + ax/2)
    x_to_minus_one_one = x -> (bx + ax - 2*x)/(-bx + ax)

    rows = 1 .+ [0;cumsum(n+1:-1:2)]
    N = (n+1)*(n+2)÷2
    pp = p(x_from_minus_one_one(prefun(Polynomial(:x₁))))
    pc = convert(ChebyshevT{FloatType, :x₁}, pp)
    c = Polynomials.coeffs(pc)
    C = zeros(eltype(c), N)
    len = min(length(c),n+1)
    C[rows[1:len]] .= c[1:len]
    C
end

function project(FloatType, p::Polynomial{T, :x₂}, n, var, prefun, boundingbox) where T
    @debug "project(poly(y)) fast path"
    if var == :x₁  # but we can see that the polynomial is in terms of x₂
        error("The function should be of :x₁")
    end
    @assert length(p) ≤ n

    ay, by = endpoints(last(boundingbox))
    y_from_minus_one_one = y -> ay - (1 + y)*(-by/2 + ay/2)
    y_to_minus_one_one = y -> (by + ay - 2*y)/(-by + ay)

    N = (n+1)*(n+2)÷2
    pp = p(y_from_minus_one_one(prefun(Polynomial(:x₂))))
    pc = convert(ChebyshevT{FloatType, :x₂}, pp)
    c = Polynomials.coeffs(pc)
    C = zeros(eltype(c), N)
    len = min(length(c),n+1)
    C[1:len] .= c[1:len]
    C
end

function project(FloatType, val::Number, n, var, prefun, boundingbox)
    @debug "project(number) fast path"
    N = (n+1)*(n+2)÷2
    C = zeros(FloatType, N)
    C[1] = val
    C
end

# Convenience function so you can do x => f instead of :x₁ => f
function subs(s::Pair{Polynomial{T, :x₁}, F}, n, nodes; prefun, boundingbox) where {T,F}
    if s[1] == Polynomials.variable(:x₁)  # and not, say, x^2 or something
        subs(:x₁=>s[2], n, nodes; prefun, boundingbox)
    else
        error("You can only substitute for x₁ or x₂.")
    end
end

# Convenience function so you can do y => f instead of :x₂ => f
function subs(s::Pair{Polynomial{T, :x₂}, F}, n, nodes; prefun, boundingbox) where {T,F}
    if s[1] == Polynomials.variable(:x₂)  # and not, say, y^2 or something
        subs(:x₂=>s[2], n, nodes; prefun, boundingbox)
    else
        error("You can only substitute for x₁ or x₂.")
    end
end

check_var_constency(var::Symbol, ::Polynomial{T, :x₁}) where T = var == :x₂
check_var_constency(var::Symbol, ::Polynomial{T, :x₂}) where T = var == :x₁
check_var_constency(var::Symbol, _) = true  # for all we know

function subs(s::Pair{Symbol, F}, n, nodes; prefun, boundingbox) where F
    if !check_var_constency(s[1], s[2])
        error("Invalid variable substitution: same variable on either side.")
    end
    if s[1] == :x₁
        subsx(s[2], n, nodes; prefun, boundingbox)
    elseif s[1] == :x₂
        subsy(s[2], n, nodes; prefun, boundingbox)
    else
        error("Symbol must be :x₁ or :x₂.")
    end
end

function subsx(funy, n, nodes; prefun, boundingbox)
    @debug "subsx(funy)"

    ax, bx = endpoints(first(boundingbox))
    x_from_minus_one_one = x -> ax - (1 + x)*(-bx/2 + ax/2)
    x_to_minus_one_one = x -> (bx + ax - 2*x)/(-bx + ax)

    ay, by = endpoints(last(boundingbox))
    y_from_minus_one_one = y -> ay - (1 + y)*(-by/2 + ay/2)
    y_to_minus_one_one = y -> (by + ay - 2*y)/(-by + ay)

    funywrap = y -> x_to_minus_one_one(funy(y_from_minus_one_one(y)))

    N = (n+1)*(n+2)÷2
    FloatType = typeof(nodes[1]*prefun(nodes[1])*funywrap(prefun(nodes[1])))
    A = zeros(FloatType, N,N)
    T(n,x) = cos(n*acos(x))
    W = reduce(hcat, T.(i,nodes) for i=0:n)
    V = reduce(hcat, T.(i,prefun.(nodes)) for i=0:n)
    Y = reduce(hcat, T.(i,funywrap.(prefun.(nodes))) for i=0:n)
    # S = reduce(hcat, Y[:,i+1] .* V[:,j+1] for i=0:n for j=0:n if i+j ≤ n)  # too slow
    idx = [(i,j) for i=0:n for j=0:n if i+j ≤ n]
    S = zeros(FloatType, length(nodes), length(idx))
    Threads.@threads for k in eachindex(idx)
        (i,j) = idx[k]
        @views S[:,k] .= Y[:,i+1] .* V[:,j+1]
    end
    A[1:n+1, :] .= 2*W'*S/length(nodes)
    A[1,:] ./= 2
    MatrixAndFunction(A, prefun, :x₂, boundingbox)
end

function subsy(funx, n, nodes; prefun, boundingbox)
    @debug "subsy(funx)"

    ax, bx = endpoints(first(boundingbox))
    x_from_minus_one_one = x -> ax - (1 + x)*(-bx/2 + ax/2)
    x_to_minus_one_one = x -> (bx + ax - 2*x)/(-bx + ax)

    ay, by = endpoints(last(boundingbox))
    y_from_minus_one_one = y -> ay - (1 + y)*(-by/2 + ay/2)
    y_to_minus_one_one = y -> (by + ay - 2*y)/(-by + ay)

    funxwrap = x -> y_to_minus_one_one(funx(x_from_minus_one_one(x)))

    N = (n+1)*(n+2)÷2
    FloatType = typeof(nodes[1]*prefun(nodes[1])*funxwrap(prefun(nodes[1])))
    A = zeros(FloatType, N,N)
    rows = 1 .+ [0;cumsum(n+1:-1:2)]
    T(n,x) = cos(n*acos(x))
    W = reduce(hcat, T.(i,nodes) for i=0:n)
    V = reduce(hcat, T.(i,prefun.(nodes)) for i=0:n)
    Y = reduce(hcat, T.(i,funxwrap.(prefun.(nodes))) for i=0:n)
    # S = reduce(hcat, V[:,i+1] .* Y[:,j+1] for i=0:n for j=0:n if i+j ≤ n) # too slow
    idx = [(i,j) for i=0:n for j=0:n if i+j ≤ n]
    S = zeros(FloatType, length(nodes), length(idx))
    Threads.@threads for k in eachindex(idx)
        (i,j) = idx[k]
        @views S[:,k] .= V[:,i+1] .* Y[:,j+1]
    end
    A[rows, :] .= 2*W'*S/length(nodes)
    A[1,:] ./= 2

    MatrixAndFunction(A, prefun, :x₁, boundingbox)
end

function subsx(val::Number, n, nodes; prefun, boundingbox)
    @debug "subsx(number) fast path"
    FloatType = eltype(nodes)
    floatval = FloatType(val)
    subsx(y->floatval, n, nodes; prefun, boundingbox)
end

function subsy(val::Number, n, nodes; prefun, boundingbox)
    @debug "subsy(number) fast path"
    FloatType = eltype(nodes)
    floatval = FloatType(val)
    subsy(x->floatval, n, nodes; prefun, boundingbox)
end

function operators(FloatType, boundingbox_, n; maxdegree::Integer=1)
    np = 2*n
    prenodes = [(2 * FloatType(k) - 1) * FloatType(π) / (2 * FloatType(np)) for k = np:-1:1]
    nodes = cos.(prenodes)

    # promote bounding box to FloatType
    boundingbox = Tuple(..(FloatType.(endpoints(i))...) for i in boundingbox_) # any tidier way?

    ax, bx = endpoints(first(boundingbox))
    scalex = (bx - ax) / 2
    shiftx = (bx + ax) / 2

    ay, by = endpoints(last(boundingbox))
    scaley = (by - ay) / 2
    shifty = (by + ay) / 2

    Dx = diffx(FloatType, n) / scalex
    Dy = diffy(FloatType, n) / scaley

    if maxdegree > 1
        Dx = tuple(Dx, (diffx(FloatType, n; degree=k) / scalex^k for k=2:maxdegree)...)
        Dy = tuple(Dy, (diffy(FloatType, n; degree=k) / scaley^k for k=2:maxdegree)...)
    end

    B = SubstitutionOperator{FloatType}(n, nodes, boundingbox)

    X = mulx(FloatType, n) * scalex + LinearAlgebra.I * shiftx
    Y = muly(FloatType, n) * scaley + LinearAlgebra.I * shifty

    Fx = f -> mulfx(FloatType, f, n, nodes, boundingbox)
    Fy = f -> mulfy(FloatType, f, n, nodes, boundingbox)

    Dx, Dy, B, X, Y, Fx, Fy
end

operators(boundingbox, n; maxdegree=1) = operators(Float64, boundingbox, n; maxdegree)
operators(n; maxdegree=1) = operators(Float64, [-1..1, -1..1], n; maxdegree)

variables() = Polynomials.variable.((:x₁, :x₂))

# TODO: think more about how best to do this type promotion stuff
getrhstype(T, x::Number, var, prefun, boundingbox) = typeof(x)
getrhstype(T, p::Polynomial, var, prefun, boundingbox) = T
function getrhstype(T, f::Function, var, prefun, boundingbox)

ax, bx = endpoints(first(boundingbox))
x_from_minus_one_one = x -> ax - (1 + x)*(-bx/2 + ax/2)

ay, by = endpoints(last(boundingbox))
y_from_minus_one_one = y -> ay - (1 + y)*(-by/2 + ay/2)

    if var == :x₁₂
        typeof(f(x_from_minus_one_one(zero(T)), y_from_minus_one_one(zero(T))))
    elseif var == :x₁
        typeof(f(x_from_minus_one_one(prefun(zero(T)))))
    else
        typeof(f(y_from_minus_one_one(prefun(zero(T)))))
    end
end

function isemptyrow(row, b)
    iszero(b) && all(iszero, row)
end

function isemptyrow(row, b, tol)
    abs(b) ≤ tol && maximum(abs, row) ≤ tol
end

# a bit hacky
stripeqns(A) = A, :x₁₂, identity, nothing
stripeqns(mfp::Pair{MatrixAndFunction{M,F}, T}) where {M,F,T} = first(mfp).A => last(mfp), first(mfp).var, first(mfp).f, first(mfp).boundingbox
# somehow we also need this next one too, in case the boundary conditions are passed by the user in a vector rather than a tuple
stripeqns(mfp::Pair{MatrixAndFunction{M}, T}) where {M,T} = first(mfp).A => last(mfp), first(mfp).var, first(mfp).f, first(mfp).boundingbox

function assemble(raweqns...; remove_empty_rows=true, pde_scaling=1)

    # the code just sort of evolved to this point, which looks a bit hacky
    eqns_and_prefuns_and_boundingboxes = stripeqns.(raweqns)
    eqns = getindex.(eqns_and_prefuns_and_boundingboxes, 1)
    vars = getindex.(eqns_and_prefuns_and_boundingboxes, 2)
    prefuns = getindex.(eqns_and_prefuns_and_boundingboxes, 3)
    boundingboxes = getindex.(eqns_and_prefuns_and_boundingboxes, 4)

    is_pde_equation = isnothing.(boundingboxes)

    ubb = unique(boundingboxes)
    # We expect to get [nothing, boundingbox] from the PDE and its boundary conditions.
    # The PDE tells you nothing about bounding box, whereas the boundary conditions
    # each smuggle in a copy of the bounding box, which should all agree.

    if length(ubb) == 1
        error("Must specify a PDE and boundary conditions.")
    end

    if length(ubb) == 2
        if !isnothing(ubb[1]) && !isnothing(ubb[2]) # one should be the PDE
            error("Internal error: confused about bounding boxes.")
        end
    end

    if length(ubb) > 2 # somehow we have more than one distinct bounding box
        error("Internal error: confused about bounding boxes.")
    end

    boundingbox = something(ubb[1], ubb[2])

    N = size(first(eqns[1]), 1)  # all matrices are the same size
    n = round(Int, sqrt(1/4 + 2*N) - 3/2)
    neq = length(eqns)

    # TODO: think more about how best to do this type promotion stuff
    MatrixFloatType = typeof(sum(zero.(eltype.(first.(eqns)))))
    RHSTypes = getrhstype.(MatrixFloatType, last.(eqns), vars, prefuns, (boundingbox,))
    RHSFloatType = typeof(zero(MatrixFloatType) + sum(zero.(RHSTypes)))

    # Assemble right hand side first, because...
    B = zeros(RHSFloatType, N, length(eqns))
    for (j,fun) in pairs(last.(eqns))
        B[:,j] = project(RHSFloatType, fun, n, vars[j], prefuns[j], boundingbox)
    end

    # ...we might want to omit rows from [A b] that are entirely zeros
    keeprow = ones(Bool, N, neq)
    # keeprow = trues(N, neq)   # BitMatrix, not thread-safe
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
    blockidxs = [0; cumsum(sum(keeprow, dims=1)[:])]
    @assert blockidxs[end] == nkeepeq
    for j = 1:neq
        Aj = first(eqns[j])
        @views A[blockidxs[j]+1:blockidxs[j+1], :] = Aj[keeprow[:,j], :]
        @views b[blockidxs[j]+1:blockidxs[j+1]] = B[keeprow[:,j], j]
        if is_pde_equation[j]
            A[blockidxs[j]+1:blockidxs[j+1], :] .*= pde_scaling
            b[blockidxs[j]+1:blockidxs[j+1]] .*= pde_scaling
        end
    end

    # Check we didn't stuff up the indexing
    # @assert isequal(reduce(vcat, M for M in first.(eqns))[keeprow[:],:], A)
    # @assert isequal(B[keeprow[:]], b)

A, b, keeprow[:], boundingbox, B[:]

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
    z = @views [F.R[1:k,1:k] \ y; zeros(Int, n-k)]   # Int so we don't promote by mistake
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

function solve(raweqns...; domain, method=:qr, cutoff=eps, remove_empty_rows=true, refine=false, scale_columns=true, pde_scaling=1)

    if method ∉ (:qr, :svd)
        error("Invalid method.  Must be :qr or :svd")
    end

    A,b,idx,boundingbox,bfull = assemble(raweqns...; remove_empty_rows, pde_scaling)

    eqns = first.(stripeqns.(raweqns))

    FloatType = eltype(A)
    rtol = cutoff == eps ? eps(FloatType)/2 : cutoff  # unit roundoff is eps/2

    if scale_columns
        colnorms = norm.(eachcol(A))
        A ./= colnorms'
    end

    F, svals = factorise!(A; method)
    c = linsolve(F, b; method, rtol)

    if scale_columns
        c ./= colnorms
    end

    CoefficientType = eltype(c)

    if refine # iterative refinement in extended precision
        # we overwrote A to save memory, but we can still get the residual from the original inputs
        bfit_extended_precision = reduce(vcat, (M*extended_precision.(c) for M in first.(eqns)))
        r = extended_precision.(b) - bfit_extended_precision[idx]
        δc = linsolve(F, r; method, rtol) # n.b. rtol is _not_ the extended precision rtol
        if scale_columns
            δc ./= colnorms
        end
        c = CoefficientType.(c + δc)
    end

    # compute fitted b based on the full inputs as a further check we didn't stuff up
    bfit = reduce(vcat, (M*c for M in first.(eqns))) # no omitting zero rows here
    # @assert all(iszero, bfit[.!idx])   # TODO: re-enable this check with a suitable tolerance

    N = sum(size(M,1) for M in first.(eqns))
    N1 = size(first(eqns[1]), 1)
    n = round(Int, sqrt(1/4 + 2*N1) - 3/2)
    maxresidual = maximum(abs, b - bfit[idx]) / maximum(abs, b) # relative to scale of b
    @info "Solve stats" n N size(A) cond(A)=svals[1]/svals[end] maxresidual
    cache = (; svals, maxresidual, residual=bfull-bfit, cutoff=svals[1]*rtol)

    reconstitute(c, boundingbox, domain; cache)
end

function reconstitute(c, boundingbox_, domain; cache)
    N = length(c)
    n = round(Int, sqrt(1/4 + 2*N) - 3/2)

    # promote bounding box to FloatType
    FloatType = eltype(c)
    boundingbox = Tuple(..(FloatType.(endpoints(i))...) for i in boundingbox_) # ugh

    PDESolution(n, c, boundingbox, domain, cache)
end

function check_condition(op, var1, var2, fun, range)
    var2 ∉ range || op(var1, fun(var2))
end

function check_condition(op, var1, var2, value::Number, range)
    var2 ∉ range || op(var1, value)
end

function raweval(u, x, y)

    ax, bx = endpoints(first(u.boundingbox))
    x_to_minus_one_one = x -> (bx + ax - 2*x)/(-bx + ax)

    ay, by = endpoints(last(u.boundingbox))
    y_to_minus_one_one = y -> (by + ay - 2*y)/(-by + ay)

    n = u.n
    dot(u.c, (cos(i*acos(x_to_minus_one_one(x)))*cos(j*acos(y_to_minus_one_one(y))) for i=0:n for j=0:n if i+j ≤ n))
end

function destructure_boundary(fun::Function, boundingbox)  # hacky way to support arbitrary mask
    fun, nothing, nothing, nothing
end

function destructure_boundary(boundary, boundingbox)
    op, varfunrange = boundary

    if varfunrange isa Tuple
        var, fun = first(varfunrange)
        range = last(varfunrange)
    else
        var, fun = varfunrange
        if var == Polynomials.variable(:x₁)
            range = last(boundingbox)
        else
            range = first(boundingbox)
        end
    end
    op, var, fun, range
end

function inrange_within_eps(x, range)
    x > range.left - eps(range.left) && x < range.right + eps(range.right)
end

function nanwrapeval(f, x, y, domain, boundingbox)

    # if x ∉ first(boundingbox) || y ∉ last(boundingbox)
    #     return NaN
    # end
    if !inrange_within_eps(x, first(boundingbox)) || !inrange_within_eps(y, last(boundingbox))
        return NaN
    end

    for boundary in domain
        op, var, fun, range = destructure_boundary(boundary, boundingbox)
        if var == Polynomials.variable(:x₁)
            if !check_condition(op, x, y, fun, range)
                return NaN
            end
        elseif var == Polynomials.variable(:x₂)
            if !check_condition(op, y, x, fun, range)
                return NaN
            end
        elseif isnothing(var) # hacky way to support arbitrary mask op(x,y)
            if !op(x, y)
                return NaN
            end
        else
            error("Boundary must be in terms of a single variable only.")
        end
    end

    f(x, y)
end

function (u::PDESolution{FloatType})(x,y; mask=true) where FloatType
    if mask
        nanwrapeval((x,y)->raweval(u, x, y), x, y, domain(u), boundingbox(u))
    else
        raweval(u, x, y)
    end
end

function __init__()
    Requires.@require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin

        function Makie.plot(
            composed_u::ComposedFunction{T,PDESolution{U,V,W}};
            divisions=256, levels=20, size=(1200, 600),
            aspect=width(first(boundingbox(composed_u.inner))) / width(last(boundingbox(composed_u.inner))),
            axis=(; aspect), diagnostics=true, mask=true,
            operator=identity, draw_boundary=false) where {T,U,V,W}

            Makie.plot(composed_u.inner; divisions, levels, size, aspect, axis, diagnostics, mask, operator, postfun=composed_u.outer)
        end

        function Makie.plot(
            u::PDESolution; divisions=256, levels=20, size=(1200, 600),
            aspect=width(first(boundingbox(u))) / width(last(boundingbox(u))),
            axis=(; aspect), diagnostics=true, mask=true,
            operator=identity, postfun=identity, draw_boundary = false
        )
            if length(divisions) == 1
                # Choose number of divisions in x and y to match the aspect ratio,
                # bounded by the passed divisions value
                if aspect > 1
                    gsx, gsy = divisions, round(Int, divisions / aspect)
                else
                    gsx, gsy = round(Int, divisions * aspect), divisions
                end
            else
                gsx, gsy = divisions
            end

            FloatType = eltype(u.c)
            bb = boundingbox(u)
            xgrid = range(first(bb).left, first(bb).right, length=gsx)
            ygrid = range(last(bb).left, last(bb).right, length=gsy)
            Z = zeros(gsx, gsy)
            Threads.@threads for i = 1:gsx
                for j = 1:gsy
                    Z[i, j] = postfun(operator(u)(xgrid[i], ygrid[j]; mask))
                end
            end
            fig = Makie.Figure(; size)
            ax, hm = Makie.heatmap(fig[1:2, 1:2], xgrid, ygrid, Z; axis)
            ax.title = "max |rₖ| / max |bₖ| = " * string(Float32(u.cache.maxresidual))
            Makie.contour!(fig[1:2, 1:2], xgrid, ygrid, Z; color=:white, levels)
            if draw_boundary
                draw_boundaries_(u.dom, u.boundingbox, max(gsx,gsy))
            end
            Makie.Colorbar(fig[1:2, 3], hm)
            if diagnostics
                N = length(u.c)
                n = round(Int, sqrt(1 / 4 + 2 * N) - 3 / 2)
                ax3, hm3 = Makie.scatter(fig[1, 4], abs.(u.c) .+ eps(FloatType), markersize=6, axis=(xlabel="𝑘", ylabel="|𝐶ₖ|", yscale=log10))
                ax3.title = "n = " * string(n) * ";  " * string(N) * " terms"
                ax2, hm2 = Makie.scatter(fig[2, 4], u.cache.svals .+ eps(FloatType), markersize=6, axis=(xlabel="𝑘", ylabel="𝑆ₖ", yscale=log10))
                if u.cache.cutoff != 0
                    Makie.lines!(fig[2, 4], [0, N], u.cache.cutoff * [1, 1]; color=:red)
                end
                ax2.title = "cond(A) ≈ " * string(Float32(u.cache.svals[1] / u.cache.svals[end]))
            end
            fig
        end

        function draw_boundaries_(dom, boundingbox, divisions) # quick and hacky
            for b = dom
                if b.x isa Tuple
                    interval = b.x[2]
                    var = b.x[1][1]
                    fun = b.x[1][2]
                else
                    var = b.x[1]
                    fun = b.x[2]
                    if var == Polynomial(:x₁)
                        interval = boundingbox[2]
                    else
                        interval = boundingbox[1]
                    end
                end
                if fun isa Number
                    bfun = t -> fun
                else
                    bfun = t -> fun(t)
                end
                t = range(interval, divisions)
                if var == Polynomial(:x₁)
                    Makie.lines!(bfun.(t), t, color=:red, linewidth=3)
                else
                    Makie.lines!(t, bfun.(t), color=:red, linewidth=3)
                end
            end
        end

    end
end

end # module EmbossPDE
