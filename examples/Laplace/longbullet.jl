using EmbossPDE
using CairoMakie

n = 60

∂x, ∂y, B = operators([-2..2, -1..1], n)
x,y = variables()

# Domain
bottom = y =>-1, -2..1
top    = y => 1, -2..1
left   = x =>-2
right  = x => y-> 1 + sqrt(1-y^2)

domain = [≥(bottom), ≤(top), ≥(left), ≤(right)]

# Equations
pde = ∂x^2 + ∂y^2 => 0
bcs = [B(bottom) => x->x, B(top) => x, B(left)*∂x => 0, B(right) => y->1+sqrt(1-y^2)]

# Solution
u = solve(pde, bcs...; domain)
@show sample = u(1,-0.25) # arbitrary point
benchmark = 1.00625295597  # see MATLAB code below
@show abs(sample-benchmark) / abs(benchmark)
plot(u)

# Lightning Laplace solver comparision https://people.maths.ox.ac.uk/trefethen/lightning.html
# u =  laplace({[1-1i 1] 1+1i -2+1i -2-1i}, {@(z) real(z), @(z) real(z), @(z) NaN*z, @(z) real(z)}, 'tol', 1e-10); u(1-0.25i)
