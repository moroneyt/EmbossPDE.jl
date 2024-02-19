using EmbossPDE
using CairoMakie

n = 60

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
bottom = y =>-1
top    = y => 1, -1/2..1/2
left  = x => (y-3)/4
right = x => (3-y)/4

domain = (≥(bottom), ≤(top), ≥(left), ≤(right))

# Equations
pde = ∂x^2 + ∂y^2 => 0
bcs = (B(bottom)*∂y => 0, B(top) => 4x^2, B(left) => y, B(right) => 1)

# Solution
u = solve(pde, bcs...; domain)
@show sample = u(0.5,-0.25) # arbitrary point
benchmark = 0.77954470798  # see MATLAB code below
@show abs(sample-benchmark) / abs(benchmark)
plot(u)

# Lightning Laplace solver benchmark https://people.maths.ox.ac.uk/trefethen/lightning.html
# u =  laplace({-1-1i, 1-1i, 0.5+1i, -0.5+1i}, {@(z) NaN*z, @(z) 0*z+1, @(z) 4*real(z).^2, @(z) imag(z)}, 'tol', 1e-10); u(0.5-0.25i)
