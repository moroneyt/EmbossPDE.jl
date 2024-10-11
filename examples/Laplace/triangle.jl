using EmbossPDE
using GLMakie

n = 40

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
bottom = y =>-1
right  = x => 1
hypot  = y => x

domain = (≥(bottom), ≤(right), ≤(hypot))

# Equations
pde = ∂x^2 + ∂y^2 => 0
bcs = (B(bottom) => x^2, B(right) => -y, B(hypot)*(∂x-∂y) => 0)

# Solution
u = solve(pde, bcs...; domain)
@show sample = u(0.5,-0.25) # arbitrary point
benchmark = 0.14436602901  # see MATLAB code below
@show abs(sample-benchmark) / abs(benchmark)
plot(u)

# Lightning Laplace solver comparision https://people.maths.ox.ac.uk/trefethen/lightning.html
# u =  laplace({1-1i, 1+1i, -1-1i}, {@(z) -imag(z), @(z) NaN*z, @(z) real(z).^2}, 'tol', 1e-10); u(0.5-0.25i)
