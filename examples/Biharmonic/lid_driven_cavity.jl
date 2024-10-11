using EmbossPDE
using GLMakie
using LinearAlgebra

n = 80

# Example 2 from Brubeck and Trefethen (2022) https://people.maths.ox.ac.uk/trefethen/stokes.pdf
∂x, ∂y, B, X, Y = operators([-1..1, -1..1], n)
x,y = variables()

# Domain
bottom = y => -1, 0..1
top    = y => 1
left   = x => -1, 0..1
right  = x => 1
curve  = y => x -> -1+sqrt(1-(x+1)^2), -1..0

domain = (≥(bottom), ≤(top), ≥(left), ≤(right), ≥(curve))

# Equations
pde = (∂x^2 + ∂y^2)^2 => 0
bcs1 = (B(bottom), B(top), B(left), B(right), B(curve)) .=> 0
bcs2 = (B(bottom)*∂y => 0, B(top)*∂y => 1, B(left)*∂x => 0, B(right)*∂x => 0, B(curve)*((X+I)*∂x + (Y+I)*∂y) => 0)

# Solution
u = solve(pde, bcs1..., bcs2...; domain)
# velsq = u -> (x,y;kwargs...) -> norm(ForwardDiff.gradient(X->u(X...;kwargs...), [x,y]))^2

exact = -0.0599323802
@show abs(u(0,0) - exact) / abs(exact)

plot(u)
# plot(u; operator=velsq, levels=0)
