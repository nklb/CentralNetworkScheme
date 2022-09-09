using CentralNetworkScheme
using Plots

f_l(u::Real) = 0.5 * u * u
f_r(u::Real) = 0.5 * u * u
λ = 1
cpl = CentralRelaxationLimit
tsteps=[ 0.2, 0.4, 0.65]#LinRange(0, 0.65, 6)# LinRange(0, 0.35, 6)
gr()

N = 200
p = Problem11(f_l, f_r, λ, cpl, Periodic)
s = constantinitialdata11(1 ,0, N)
s.u₁= .5 .+ .5 * sin.(pi .* (1 .+ collect(s.xc[1:N-1])))
s.u₂= .5 .+ .5 * sin.(pi .* (1 .+ collect(s.xc[N:end])))

y = ones(3)
println("Computing solution for ", cpl)
title = Plots.scatter(y, marker=0, markeralpha=0,
                      axis=false, grid=false, leg=false)
plots = []
for t in tsteps
    global s = centralcoupling(p, s, t, FirstOrder, .9)
    mass = dx(s) * sum(u(s))
    push!(plots, scatter(s.xc, u(s), markershape=:plus, label="",
                         title=string("t=", t)))
end

s = constantinitialdata11(1 ,0, N)
s.u₁= .5 .+ .5 * sin.(pi .* (1 .+ collect(s.xc[1:N-1])))
s.u₂= .5 .+ .5 * sin.(pi .* (1 .+ collect(s.xc[N:end])))
k = 0
for t in tsteps
    global s = centralcoupling(p, s, t, SecondOrder, .2)
    mass = dx(s) * sum(u(s))
    push!(plots, scatter(s.xc, u(s), markershape=:plus, label=""))
end

Plt = plot(title, plot(plots ..., layout=(2,3), size=(1800, 900)),
           layout=grid(2,1,heights=[0.1,0.9]))
display(Plt)
