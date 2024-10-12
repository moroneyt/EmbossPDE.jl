using EmbossPDE
using GLMakie
using ForwardDiff: derivative

n = 60

∂x, ∂y, B, X, Y, Fx, Fy = operators([0..1, 0..1], n)
x,y = variables()

f  = x -> x + 0.05*sinpi(4x)  # but not too wavy
f′ = x -> derivative(f,x)

# Domain
bottom = y => 0
right  = x => 1
hypot  = y => f

domain = (≥(bottom), ≤(right), ≤(hypot))

# Equations
pde = ∂x^2 + ∂y^2 => 0
bcs = (B(bottom) => x, B(right) => 1-y, B(hypot)*(Fx(f′)*∂x-∂y) => 0)

# Solution
u = solve(pde, bcs...; domain)
plot(u)