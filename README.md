# EmbossPDE

EMbedded BOundary Spectral Solver for PDEs in two dimensions.

⚠️This is a work in progress.

## Example
Laplace's equation on a triangle with Dirichlet and Neumann boundary conditions.
```
# Only needed the first time
using Pkg
pkg"activate --temp"
pkg"add https://github.com/moroneyt/EmbossPDE.jl, GLMakie"
using EmbossPDE
using GLMakie: plot

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
plot(u)

```
![Triangle solution](https://github.com/moroneyt/EmbossPDE.jl/blob/main/triangle.png)

### Notes
* See `examples` directory for more

[![Build Status](https://github.com/moroneyt/EmbossPDE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/moroneyt/EmbossPDE.jl/actions/workflows/CI.yml?query=branch%3Amain)
