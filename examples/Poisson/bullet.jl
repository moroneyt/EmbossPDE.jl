using EmbossPDE
using CairoMakie

n = 40

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
bottom = y =>-1, -1..0
top    = y => 1, -1..0
left   = x =>-1
curve  = x => y->sqrt(1-y^2)

domain = (≥(bottom), ≤(top), ≥(left), ≤(curve))

# Equations
pde = ∂x^2 + ∂y^2 => 1
bcs = (B(bottom), B(top), B(left), B(curve)) .=> 0

# Solution
u = solve(pde, bcs...; domain)
plot(u)