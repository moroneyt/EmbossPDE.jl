using EmbossPDE
using CairoMakie

n = 40

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
C = sqrt(1-(-1+0.5)^2)
bottom = y => -0.5, -C..C
left   = x => y->-sqrt(1-y^2), -0.5..1
right  = x => y-> sqrt(1-y^2), -0.5..1

domain = (≥(y=>-0.5), ≥(left), ≤(right))

# Equations
pde = ∂x^2 + ∂y^2 => 1
bcs = (B(bottom), B(left), B(right)) .=> 0

# Solution
u = solve(pde, bcs...; domain)
plot(u)