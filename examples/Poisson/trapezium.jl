using EmbossPDE
using GLMakie

n = 40

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
bottom = y =>-1
top    = y => 1, -1/2..1/2
left   = x => (y-3)/4
right  = x => (3-y)/4

domain = (≥(bottom), ≤(top), ≥(left), ≤(right))

# Equations
pde = ∂x^2 + ∂y^2 => 1
bcs = (B(bottom), B(top), B(left), B(right)) .=> 0

# Solution
u = solve(pde, bcs...; domain)
plot(u)
