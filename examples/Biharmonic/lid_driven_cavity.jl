using EmbossPDE
using GLMakie

n = 80

∂x, ∂y, B = operators([-1..1, -1..1], n)
x,y = variables()

# Domain
bottom = y => -1
top    = y => 1
left   = x => -1
right  = x => 1

domain = (≥(bottom), ≤(top), ≥(left), ≤(right))

# Equations
pde = (∂x^2 + ∂y^2)^2 => 0
bcs1 = (B(bottom)=> 0, B(top)=> 0, B(left)=> 0, B(right)=> 0)
bcs2 = (B(bottom)*∂y => 0, B(top)*∂y => -1, B(left)*∂x => 0, B(right)*∂x => 0)

# Solution
u = solve(pde, bcs1..., bcs2...; domain)

exact = 0.11790231118443;
@show abs(u(0,0) - exact) / exact

plot(u)

