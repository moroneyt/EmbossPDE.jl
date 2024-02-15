using EmbossPDE
using CairoMakie

n = 40

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
top    = y => 1
bottom = y =>-1, -1..0
left   = x =>-1
right  = x => 1,  0..1
edge   = y => x-1,  0..1

domain = (≥(bottom), ≤(top), ≥(left), ≤(right), ≥(edge))

# Equations
pde = ∂x^2 + ∂y^2 => 1
bcs = (B(bottom), B(top), B(left), B(right), B(edge)) .=> 0

# Solution
u = solve(pde, bcs...; domain)
plot(u)