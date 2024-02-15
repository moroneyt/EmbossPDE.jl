using EmbossPDE
using CairoMakie

n = 40

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
bottom = y =>-1
right  = x => 1
hypot  = y => 1/2*(x+1)^2-1

domain = (≥(bottom), ≤(right), ≤(hypot))

# Equations
pde = ∂x^2 + ∂y^2 => 1
bcs = (B(bottom), B(right), B(hypot)) .=> 0

# Solution
u = solve(pde, bcs...; domain)
plot(u)