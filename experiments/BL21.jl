using CentralNetworkScheme
using Plots
 
a = .5
f(u::Real) =  u * u / (u * u + a * (1 - u) * (1 - u)) 

λ = 2.5
tsteps= [0, 0.25, 0.45, 0.8] #LinRange(0, 0.7, 4)
bc = HomNeumann
gr()

cpl = CentralRelaxationLimit

p = Problem21(f, f, f, λ, cpl, bc)
s = constantinitialdata21(0, .16, 0, 300) 
s.u₁[collect(s.xc₁) .< -.5] .= 1


plots = []

for t in tsteps
    global s = centralcoupling(p, s, t, SecondOrder, 0.49)
    incoming = plot(collect(s.xc₁), s.u₁, color = "red", label = "",
                    title = string("t=", t, " incoming edges"), ylims = (0, 1.1))
    plot!(incoming, s.xc₂, s.u₂, color = "blue", linestyle = :dash, label = "")
    push!(plots, incoming)
    push!(plots, plot(s.xc₃, s.u₃, color = "black", label = "",
                      title = string("outgoing edge"), ylims = (0, 1.1)))
end

Plt = plot(plots..., layout = (4, 2), size = (1800, 1500))
display(Plt)
