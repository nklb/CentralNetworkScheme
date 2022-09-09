using CentralNetworkScheme
using Plots
 
f(u::Real) = u * (1-u)
g(u::Real) = u * (1-u/1.2)

λ = 1
tsteps= [0, .25, .75, 2]
bc = Outgoing
gr()

cpl = CentralRelaxationLimit

β = .2
p = Problem21(f, f, g, λ, cpl, bc, β, true)
s = constantinitialdata21(0.07, 0.15, 0.2)

y = ones(3)
plots = []

for t in tsteps
    global s = centralcoupling(p, s, t, FirstOrder, 0.49)
    mass = (s.xc₁[2] - s.xc₁[1]) * sum(s.u₁) + (s.xc₂[2] - s.xc₂[1]) * sum(s.u₂) +
        (s.xc₃[2] - s.xc₃[1]) * sum(s.u₃)
    incoming = plot(collect(s.xc₁), s.u₁, color = "red", label = "",
                    title = string("t=", t, " incoming roads"), ylims = (0, 0.4))
    plot!(incoming, s.xc₂, s.u₂, color = "blue", linestyle = :dash, label = "")
    push!(plots, incoming)
    push!(plots, plot(s.xc₃, s.u₃, color = "black", label = "",
                      title = string("outgoing road"), ylims = (0, 0.4)))
end

title = Plots.scatter(y, marker = 0, markeralpha = 0,
                                 axis = false, grid = false, leg = false)
Plt = plot(title, plot(plots..., layout = (4, 2), size = (1800, 1500)),
           layout = grid(2, 1, heights = [0.1, 0.9]))
display(Plt)           

