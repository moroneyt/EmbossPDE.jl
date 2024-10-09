using EmbossPDE
using GLMakie

n = 60

∂x, ∂y, B, X, Y = operators([-1..1, 0..1], n)
x,y = variables()

# Domain
top    = y => x-> sqrt(1-x^2)
bottom = y => 0

domain = (≥(bottom), ≤(top))

# Equations
pde = ∂x^2 + ∂y^2 => 1
bcs = (B(bottom) => x, B(top)*(X*∂x+Y*∂y) => 0)

## Solution
u = solve(pde, bcs...; domain)
@show u(0.25,0.25√3)  # compare with polar.jl
plot(u)