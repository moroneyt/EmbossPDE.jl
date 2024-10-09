using EmbossPDE
using GLMakie

n = 80

M = 0.8
MM = 1 - M

∂x, ∂y, B = operators([0..5, -1..1], n)
x,y = variables()

# Domain
bottom = y => x->-M + MM/2*sinpi(x) + MM/4*sinpi(2*x)
top    = y => x->M - MM/3*sinpi(x) + MM/5*sinpi(2*x)
left   = x => 0, -M..M
right  = x => 5, -M..M

domain = [≥(bottom), ≤(top), ≥(left), ≤(right)]

pde = ∂x^2 + ∂y^2 => 1

bcs = (B(bottom) => 0, B(top) => 0, B(left) - B(right) => 0, B(left)*∂x - B(right)*∂x => 0)

# Solution
u = solve(pde, bcs...; domain)
plot(u, divisions=512)

