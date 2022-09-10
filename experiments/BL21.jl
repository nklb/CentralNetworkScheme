using CentralNetworkScheme
using Plots

# setup of flux functions on all 3 edges
a = .5
f(u::Real) =  u * u / (u * u + a * (1 - u) * (1 - u)) 

# relaxation speed and time instances for plotting
λ = 2.5
tsteps = [0, 0.25, 0.45, 0.8]

# boundary condition, choose between ZeroFlux HomNeumann Outgoing Periodic
bc = HomNeumann 

# setup of the problem and the initial condition
p = Problem21(f, f, f, λ, CentralRelaxationLimit, bc)
s = constantinitialdata21(0, .16, 0, 300) 
s.u₁[collect(s.xc₁) .< -.5] .= 1

# initialize subplots array and gr plotting framework
plots = []
gr()

# simulation and plotting loop
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
