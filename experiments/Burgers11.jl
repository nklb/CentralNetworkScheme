using CentralNetworkScheme
using Plots

# setup of flux functions on incoming and outgoing edge
f_l(u::Real) = 0.5 * u * u
f_r(u::Real) = 0.5 * u * u

# relaxation speed, time instances for plotting and mesh cells/edge
λ = 1
tsteps=[ 0.2, 0.4, 0.65]
N = 200

# setup of the problem and the initial condition
p = Problem11(f_l, f_r, λ, CentralRelaxationLimit, Periodic)
s = constantinitialdata11(1 ,0, N)
s.u₁= .5 .+ .5 * sin.(pi .* (1 .+ collect(s.xc[1:N-1])))
s.u₂= .5 .+ .5 * sin.(pi .* (1 .+ collect(s.xc[N:end])))

# setup subplots array and gr plotting framework
plots = []
gr()

# simulation and plotting loop for the first order scheme
for t in tsteps
    global s = centralcoupling(p, s, t, FirstOrder, .9)
    push!(plots, scatter(s.xc, u(s), markershape=:plus, label="",
                         title=string("t=", t)))
end

# reset problem and initial condition
s = constantinitialdata11(1 ,0, N)
s.u₁= .5 .+ .5 * sin.(pi .* (1 .+ collect(s.xc[1:N-1])))
s.u₂= .5 .+ .5 * sin.(pi .* (1 .+ collect(s.xc[N:end])))

# simulation and plotting loop for the second order scheme
for t in tsteps
    global s = centralcoupling(p, s, t, SecondOrder, .2)
    push!(plots, scatter(s.xc, u(s), markershape=:plus, label=""))
end

Plt = plot(plots ..., layout=(2,3), size=(1800, 900))
display(Plt)
