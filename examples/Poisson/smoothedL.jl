using EmbossPDE
using CairoMakie

n = 40

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
C = 0.5
bottom = y =>-1
top    = y => 1, -1..0
left   = x =>-1
right  = x => 1, -1..0
right2 = x => 0, C .. 1
top2   = y => 0, C .. 1
curve  = y => x->C-sqrt(C^2-(x-C)^2), 0..C

domain = (≥(bottom), ≤(top), ≥(left), ≤(right), ≤(right2), ≤(top2), ≤(curve))

# Equations
pde = ∂x^2 + ∂y^2 => 1
bcs = (B(bottom), B(top), B(left), B(right), B(right2), B(top2), B(curve)) .=> 0

# Solution
u = solve(pde, bcs...; domain)
plot(u)