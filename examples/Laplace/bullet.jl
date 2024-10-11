using EmbossPDE
using GLMakie

n = 60

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
bottom = y =>-1, -1..0
top    = y => 1, -1..0
left   = x =>-1
right  = x => y->sqrt(1-y^2)

domain = (≥(bottom), ≤(top), ≥(left), ≤(right))

# Equations
pde = ∂x^2 + ∂y^2 => 0
bcs = (B(bottom) => x, B(top) => x, B(left)*∂x => 0, B(right) => y->sqrt(1-y^2))

# Solution
u = solve(pde, bcs...; domain)
@show sample = u(0.5,-0.25) # arbitrary point
benchmark = 0.54953078508  # see MATLAB code below
@show abs(sample-benchmark) / abs(benchmark)
plot(u)

# Lightning Laplace solver benchmark https://people.maths.ox.ac.uk/trefethen/lightning.html
# u =  laplace({[0-1i 1] 0+1i -1+1i -1-1i}, {@(z) real(z), @(z) real(z), @(z) NaN*z, @(z) real(z)}, 'tol', 1e-10); u(0.5-0.25i)
