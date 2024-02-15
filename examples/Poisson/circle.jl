using EmbossPDE
using CairoMakie

n = 40

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
top    = y => x-> sqrt(1-x^2)
bottom = y => x->-sqrt(1-x^2)

domain = (≥(bottom), ≤(top))

# Equations
pde = ∂x^2 + ∂y^2 => 1
bcs = (B(bottom), B(top)) .=> 0

# Solution
u = solve(pde, bcs...; domain)
plot(u)