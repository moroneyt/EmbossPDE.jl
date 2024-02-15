using EmbossPDE
using CairoMakie

n = 40

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
bottom = y =>-1
top    = y => 1
left   = x =>-1
right  = x => 1

domain = (≥(bottom), ≤(top), ≥(left), ≤(right))

# Equations
pde = ∂x^2 + ∂y^2 => 1
bcs = (B(bottom), B(top), B(left), B(right)) .=> 0

# Solution
u = solve(pde, bcs...; domain)
plot(u)