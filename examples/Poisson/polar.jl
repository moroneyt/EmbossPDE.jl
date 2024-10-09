using EmbossPDE
using GLMakie

n = 80

∂r, ∂θ, B, R, Θ = operators([0..1, 0..π], n)
r,θ = variables()

# Domain
left   = r => 0
right  = r => 1
bottom = θ => 0
top    = θ => π

domain = (≥(bottom), ≤(top), ≥(left), ≤(right))

# Equations
pde = R*∂r*(R*∂r) + ∂θ^2 => r^2
bcs = (B(left)*∂r => 0, B(right)*∂r => 0, B(bottom) => r, B(top) => -r)

# Solution
u = solve(pde, bcs...; domain)
@show u(0.5,π/3)  # compare with semicircle.jl
plot(u, aspect=2)