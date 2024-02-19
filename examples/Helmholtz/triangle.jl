using EmbossPDE
using CairoMakie
using LinearAlgebra

n = 60

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
bottom = y =>-1
right  = x => 1
hypot  = y => x

domain = (≥(bottom), ≤(right), ≤(hypot))

# Equations
pde = ∂x^2 + ∂y^2 + 50^2/4*I => 1
bcs = (B(bottom), B(right), B(hypot)) .=> 0

# Solution
u = solve(pde, bcs...; domain)
plot(u; levels=0)