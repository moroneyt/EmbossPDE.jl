using EmbossPDE
using GLMakie
using LinearAlgebra
using ForwardDiff: gradient, hessian
shg = fig -> display(GLMakie.Screen(), fig)

n = 60  # polynomial degree

∂x, ∂y, B = operators([-1..1, -1..1], n)
x,y = variables()

# Domain
bottom = y =>-1
right  = x => 1
hypot  = y => x

domain = (≥(bottom), ≤(right), ≤(hypot))

# Equations
pde = ∂x^2 + ∂y^2 => 0
bc1 = B(bottom) => x^2
bc2 = B(right) => -y
bc3 = B(hypot)*(∂x-∂y) => 0

# Solution
u = solve(pde, bc1, bc2, bc3; domain)
plot(u) |> shg

## Residuals
∇(u)  = (x,y) -> gradient(X->u(X...), [x,y])
∇²(u) = (x,y) -> tr(hessian(X->u(X...), [x,y]))

# PDE
plot(∇² ∘ u; diagnostics=false, axis=(title="∇²u",)) |> shg

# Boundaries
fig = Figure()
lines(fig[1, 1], -1..1, x -> u(x,-1) - x^2,    axis=(title="bottom",))
lines(fig[2, 1], -1..1, y -> u(1,y) + y,       axis=(title="right",))
lines(fig[3, 1], -1..1, a -> ∇(u)(a,a)⋅[1,-1], axis=(title="hypot",))
fig |> shg
