using EmbossPDE
using GLMakie

n = 60

∂x, ∂y, B = operators([-1..1, -1..1], n)
x,y = variables()

# Domain
north     = y => (x+1/2)*(x-1/2)+1,        -1/2..1/2
northeast = x => -y^2/2+1,                    0..1
southeast = x => y->sinpi(y/2)/2+1,          -1..0
south     = y => -1,                       -1/2..1/2
southwest = x => y->-0.1*sinpi(2*y)-y/2-1,   -1..0
northwest = x => y->sqrt(y)/2-1,              0..1

domain = [.≥([south, southwest, northwest]); .≤([north, northeast, southeast])]

# Equations
pde = ∂x^2 + ∂y^2 => 1
bcs = B.([south, north, southwest, northwest, southeast, northeast]) .=> 0

## Solution
u = solve(pde, bcs...; domain)
plot(u; divisions=512)