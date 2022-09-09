using CentralNetworkScheme
using Plots
 
a = .5
f(u::Real) =  u * u / (u * u + a * (1 - u) * (1 - u)) 
g(u::Real) =  u * u / (u * u + a * (1 - u) * (1 - u)) 

λ = 2.5
tsteps= [0, 0.25, 0.45, 0.8] #LinRange(0, 0.7, 4)
bc = HomNeumann
gr()

cpl = CentralRelaxationLimit

 p = Problem21(f, f, g, λ, cpl, bc)
 s = constantinitialdata21(0, .16, 0, 300) 
s.u₁[collect(s.xc₁) .< -.5] .= 1

y = ones(3)

plots = []

for t in tsteps
    global s = centralcoupling(p, s, t, SecondOrder, 0.49)
    mass = (s.xc₁[2] - s.xc₁[1]) * sum(s.u₁) + (s.xc₂[2] - s.xc₂[1]) * sum(s.u₂) +
        (s.xc₃[2] - s.xc₃[1]) * sum(s.u₃)
    incoming = plot(collect(s.xc₁), s.u₁, color = "red", label = "",
                    title = string("t=", t, " incoming edges"), ylims = (0, 1.1))
    plot!(incoming, s.xc₂, s.u₂, color = "blue", linestyle = :dash, label = "")
    push!(plots, incoming)
    push!(plots, plot(s.xc₃, s.u₃, color = "black", label = "",
                      title = string("outgoing edge"), ylims = (0, 1.1)))
end
 title = Plots.scatter(y, marker = 0, markeralpha = 0,
                             axis = false, grid = false, leg = false)
Plt = plot(title, plot(plots..., layout = (4, 2), size = (1800, 1500)),
           layout = grid(2, 1, heights = [0.1, 0.9]))
display(Plt)
