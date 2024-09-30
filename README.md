# EmbossPDE

EMbedded BOundary Spectral Solver for PDEs in two dimensions.

⚠️This is a work in progress.

## Example
Laplace's equation on a triangle with Dirichlet and Neumann boundary conditions.
```
using Pkg
pkg"activate --temp"
pkg"add https://github.com/moroneyt/EmbossPDE.jl, GLMakie"
using EmbossPDE
using GLMakie: plot

n = 30  # polynomial degree

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
bottom = y =>-1
right  = x => 1
hypot  = y => x

domain = (≥(bottom), ≤(right), ≤(hypot))

# Equations
pde = ∂x^2 + ∂y^2 => 0
bcs = (B(bottom) => x^2, B(right) => -y, B(hypot)*(∂x-∂y) => 0)

# Solution
u = solve(pde, bcs...; domain)
plot(u)
```

### Notes
* Domain is a simply-connected subset of $[-1,1]^2$
* Linear problems only
* See `examples` directory for more

[![Build Status](https://github.com/moroneyt/EmbossPDE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/moroneyt/EmbossPDE.jl/actions/workflows/CI.yml?query=branch%3Amain)
