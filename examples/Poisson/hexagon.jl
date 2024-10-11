using EmbossPDE
using GLMakie

n = 40

∂x, ∂y, B = operators(n)
x,y = variables()

# Domain
south = y =>-1, -1/2..1/2
north = y => 1, -1/2..1/2
southwest = x =>-y/2-1, -1..0
northwest = x => y/2-1,  0..1
southeast = x => y/2+1, -1..0
northeast = x =>-y/2+1,  0..1
domain = (≥(south), ≤(north), ≥(southwest), ≥(northwest), ≤(southeast), ≤(northeast))

# Equations
pde = ∂x^2 + ∂y^2 => 1
bcs = (B(south), B(north), B(southwest), B(northwest), B(southeast), B(northeast)) .=> 0

# Solution
u = solve(pde, bcs...; domain)
plot(u)
