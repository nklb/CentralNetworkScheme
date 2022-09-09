module CentralNetworkScheme

using LinearAlgebra

export Coupling, BoundaryCondition, Scheme, NumSolution11, NumSolution21, Problem11, Problem21
export constantinitialdata11, constantinitialdata21, centralcoupling, centraluncoupled
export u, dx

@enum Coupling CentralRelaxationLimit LocalRelaxation GodlewskiRaviart TrafficFlowMaximization
@enum BoundaryCondition ZeroFlux HomNeumann Outgoing Periodic
@enum Scheme FirstOrder SecondOrder


mutable struct NumSolution11
    u₁::Vector{<:Real}
    u₂::Vector{<:Real}
    N::Integer
    t::Real
    xc::LinRange
end


mutable struct NumSolution21
    u₁::Vector{<:Real}
    u₂::Vector{<:Real}
    u₃::Vector{<:Real}
    N::Integer
    t::Real
    xc₁::LinRange
    xc₂::LinRange
    xc₃::LinRange
end


struct Problem11
    f₁
    f₂
    λ::Real
    cpl::Coupling
    bc::BoundaryCondition
    fixeddt::Real
    nodeslope::Bool
end


mutable struct Problem21
    f₁
    f₂
    f₃
    λ::Real
    cpl::Coupling
    bc::BoundaryCondition
    β::Real # right of way parameter
    nodeslope::Bool
end


function u(s::NumSolution11)
    [s.u₁; s.u₂]
end


function dx(s::NumSolution11)
    dx = (s.xc.stop - s.xc.start) / s.xc.lendiv
end


function constantinitialdata11(u_left::Real, u_right::Real, N::Integer=201)
    u₁ = zeros(N-1)
    u₂ = zeros(N-1)
    u₁ .= u_left
    u₂ .= u_right
    t = 0
    xi = LinRange(-1, 1, 2 * N - 1)
    xc = .5 * (xi[2:end] + xi[1:end-1])
    NumSolution11(u₁, u₂, N, t, xc)
end


function constantinitialdata21(u_1::Real, u_2::Real, u_3::Real, N::Integer=201)
    u₁ = zeros(N-1)
    u₂ = zeros(N-1)
    u₃ = zeros(N-1)
    u₁ .= u_1
    u₂ .= u_2
    u₃ .= u_3
    t = 0
    xc₁ = LinRange(-1, 0, N-1)
    xc₂ = LinRange(-1, 0, N-1)
    xc₃ = LinRange(0, 1, N-1)
    NumSolution21(u₁, u₂, u₃, N, t, xc₁, xc₂, xc₃)
end


function Problem11(f₁, f₂, λ::Real, cpl::Coupling, bc::BoundaryCondition)
    Problem11(f₁, f₂, λ, cpl, bc, 0, false)
end


function numflux(u_l::Real, u_r::Real, v_l::Real, v_r::Real, s_l⁺::Real, s_r⁻::Real,
                 λ::Real)
    0.5 * (v_r + v_l) - 0.5 * λ * (u_r - u_l) - .5 * (s_r⁻ - s_l⁺)
end


function coupling(p::Problem11, u_trL::Real, u_trR::Real)
    if p.cpl == CentralRelaxationLimit
        u_R = (u_trL + u_trR) / 2 + (p.f₁(u_trL) - p.f₂(u_trR)) / (2 * p.λ)
        u_L = u_R
        v_R = (p.f₁(u_trL) + p.f₂(u_trR)) / 2 + (p.λ / 2) * (u_trL - u_trR)
        v_L = v_R
    elseif p.cpl == LocalRelaxation
        u_R = (u_trL + u_trR) / 2 + (p.f₁(u_trL) - p.f₂(u_trR)) / (2 * p.λ)
        u_L = u_R
        v_R = p.f₁(u_R)
        v_L = p.f₂(u_L)
    elseif p.cpl == GodlewskiRaviart
        u_R = u_trR
        u_L = u_trL
        v_R = p.f₁(u_R)
        v_L = p.f₂(u_L)
    end
    u_R, u_L, v_R, v_L
end


function coupling(p::Problem21, u_tr1::Real, u_tr2::Real, u_tr3::Real)
    if p.cpl == CentralRelaxationLimit

        S = Diagonal([-1, -1, 1])
        u = [u_tr1, u_tr2, u_tr3]
        v = [p.f₁(u_tr1), p.f₂(u_tr2), p.f₃(u_tr3)]
        β = [v[2], -v[1], 0]
        b = [v[3] - v[1] - v[2], p.λ^2 * (u[3] - u[1] - u[2]), -(β[1] * v[1] + β[2] * v[2] + β[3] * v[3])]

        M = p.λ * [1 1 -1; -p.λ -p.λ -p.λ; transpose(β)]
        Mₛ = M * S

        uc = Mₛ\(Mₛ * u + b)
        vc = M\(M * v + p.λ*b)

        u_R1 = uc[1]
        u_R2 = uc[2]
        u_L3 = uc[3]
        v_R1 = vc[1]
        v_R2 = vc[2]
        v_L3 = vc[3]

    elseif p.cpl == TrafficFlowMaximization
        # if this is chosen, it's assumed that flux functions are f₁(ρ)=f₂(ρ)=ρ(1-ρ)
        # and f₃(ρ)=ρ(1-ρ/1.2) 
        r₁=1
        r₂=1.2
        demand(ρ) = ρ <= r₁/2 ? ρ*(1-ρ/r₁) : r₁/2
        supply(ρ) = ρ <= r₂/2 ? r₂/2 : ρ*(1-ρ/r₂)

        d1 = demand(u_tr1)
        d2 = demand(u_tr2)
        maxin = d1 + d2 
        maxout = supply(u_tr3)
        v_L3 = minimum([maxin, maxout])
        println("total demand: ", maxin,(d1, d2), ", total supply: ", maxout)
        if maxin <= maxout # free flow case
            v_R1 = demand(u_tr1)
            v_R2 = demand(u_tr2)
        else
            v_R1 = p.β * maxout
            v_R2 = (1 - p.β) * maxout  
            if v_R1 > d1
                v_R1 = d1
                v_R2 = maxout - d1
            elseif v_R2 > d2
                v_R2 = d2
                v_R1 = maxout - d2
            end
        end
        u_R1 = 0
        u_R2 = 0
        u_L3 = 0
    end
    u_R1, u_R2, u_L3, v_R1, v_R2, v_L3
end


function mcslopes(v::Vector{<:Real})
    dv = 2 * diff(v)
    central = (v[3:end] - v[1:end-2]) / 2
    sdv = sign.(dv)

    (sdv[2:end] .== sdv[1:end-1]) .* sdv[2:end] .*
        min.(abs.(dv[1:end-1]), abs.(central), abs.(dv[2:end]))
end


function charslopes(u::Vector{<:Real}, v::Vector{<:Real}, λ::Real)
    mcslopes(.5 * (v - λ .* u)), mcslopes(.5 * (v + λ .* u))
end


function ghostdata(s::NumSolution11, p::Problem11)
    if p.bc == Periodic
        ughost = (s.u₂[end], s.u₁[1])
        vghost = (p.f₂(s.u₂[end]), p.f₁(s.u₁[1]))
    else
        ughost = (s.u₁[1], s.u₂[end])
        vghost = (p.f₁(s.u₁[1]), p.f₂(s.u₂[end]))
    end
    ughost, vghost
end

    
function centralcoupling(p::Problem11, s::NumSolution11, T::Real, sc::Scheme=FirstOrder, CFL::Real = 0.9)
    dx = (s.xc.stop - s.xc.start) / s.xc.lendiv
    if p.fixeddt > 0
        dt = p.fixeddt
    else
        dt = CFL * dx / p.λ
    end
    N = s.N
    warned = false
    nfx₁ = zeros(N)
    nfx₂ = zeros(N)

    s₁⁻ = zeros(N-1)
    s₁⁺ = zeros(N-1)
    s₂⁻ = zeros(N-1)
    s₂⁺ = zeros(N-1)
    sghost = (0, 0)

    while s.t<T

        u_R, u_L, v_R, v_L = coupling(p, s.u₁[end], s.u₂[1])

        ughost, vghost = ghostdata(s, p)
            
        if sc == SecondOrder
            s₁⁻, s₁⁺ = charslopes(vcat(ughost[1], s.u₁, u_R), vcat(vghost[1], p.f₁.(s.u₁), v_R), p.λ)
            s₂⁻, s₂⁺ = charslopes(vcat(u_L, s.u₂, ughost[2]), vcat(v_L, p.f₂.(s.u₂), vghost[2]), p.λ)
            if p.bc == Periodic
                sghost = (s₂⁺[end], s₁⁻[1])
            end
        end

        if p.bc in [HomNeumann, Periodic]
            nfx₁[1] = numflux(ughost[1], s.u₁[1], vghost[1], p.f₁(s.u₁[1]), sghost[1], s₁⁻[1], p.λ)
            nfx₂[N] = numflux(s.u₂[end], ughost[2], p.f₂(s.u₂[end]), vghost[2], s₂⁺[end], sghost[2], p.λ)
        end

        if !p.nodeslope
             s₁⁻[end]= 0
             s₂⁺[1] = 0
        end
        
        # fluxes left and right
        for k=2:N-1
            nfx₁[k] = numflux(s.u₁[k-1], s.u₁[k], p.f₁(s.u₁[k-1]), p.f₁(s.u₁[k]), s₁⁺[k-1], s₁⁻[k], p.λ)
            nfx₂[k] = numflux(s.u₂[k-1], s.u₂[k], p.f₂(s.u₂[k-1]), p.f₂(s.u₂[k]), s₂⁺[k-1], s₂⁻[k], p.λ)
        end

        # @show s₁⁺[end], s₂⁻[1]
        # @show s₁⁻[end], s₂⁺[1]
        # fluxes center
        nfx₁[N] = numflux(s.u₁[end], u_R, p.f₁(s.u₁[end]), v_R, s₁⁺[end], 0, p.λ)
        nfx₂[1] = numflux(u_L, s.u₂[1], v_L, p.f₂(s.u₂[1]), 0, s₂⁻[1], p.λ)
        
        if abs(nfx₁[N] - nfx₂[1]) > 1e-12 && !warned
            println("Scheme not conservative")
            warned = true
        end
        
        # state update
        s.u₁ .-= dt/dx * diff(nfx₁)
        s.u₂ .-= dt/dx * diff(nfx₂)
        s.t += dt
        println("t=", s.t, ", TV=", sum(abs.(diff(u(s)))))
    end
    s
end


function centralcoupling(p::Problem21, s::NumSolution21, T::Real, sc::Scheme=FirstOrder, CFL::Real = 0.9)
    dx = (s.xc₁.stop - s.xc₁.start) / s.xc₁.lendiv
    dt = CFL * dx / p.λ
    N = s.N
    warned = false
    nfx₁ = zeros(N)
    nfx₂ = zeros(N)
    nfx₃ = zeros(N)

    s₁⁻ = zeros(N-1)
    s₁⁺ = zeros(N-1)
    s₂⁻ = zeros(N-1)
    s₂⁺ = zeros(N-1)
    s₃⁻ = zeros(N-1)
    s₃⁺ = zeros(N-1)

    while s.t<T
        
        u_R1, u_R2, u_L3, v_R1, v_R2, v_L3 = coupling(p, s.u₁[end], s.u₂[end], s.u₃[1])

        if sc == SecondOrder
            s₁⁻, s₁⁺ = charslopes(vcat(s.u₁[1], s.u₁, u_R1), vcat(p.f₁(s.u₁[1]), p.f₁.(s.u₁), v_R1), p.λ)
            s₂⁻, s₂⁺ = charslopes(vcat(s.u₂[1], s.u₂, u_R2), vcat(p.f₂(s.u₂[1]), p.f₂.(s.u₂), v_R2), p.λ)
            s₃⁻, s₃⁺ = charslopes(vcat(u_L3, s.u₃, s.u₃[end]), vcat(v_L3, p.f₃.(s.u₃), p.f₃(s.u₃[end])), p.λ)
        end
        
        #@show s₂⁺[end], s₁⁺[end], s₃⁻[1]
        if p.bc == HomNeumann
            nfx₁[1] = numflux(s.u₁[1], s.u₁[1], p.f₁(s.u₁[1]), p.f₁(s.u₁[1]), 0, s₁⁻[1], p.λ)
            nfx₂[1] = numflux(s.u₂[1], s.u₂[1], p.f₂(s.u₂[1]), p.f₂(s.u₂[1]), 0, s₂⁻[1], p.λ)
        end
        if p.bc in [HomNeumann, Outgoing]
            nfx₃[N] = numflux(s.u₃[end], s.u₃[end], p.f₃(s.u₃[end]), p.f₃(s.u₃[end]), s₃⁺[end], 0, p.λ)
        end

        if !p.nodeslope
            s₁⁻[end]= 0
            s₂⁻[end]= 0
            s₃⁺[1] = 0
        end

        # fluxes away from the junction
        for k=2:N-1
            nfx₁[k] = numflux(s.u₁[k-1], s.u₁[k], p.f₁(s.u₁[k-1]), p.f₁(s.u₁[k]), s₁⁺[k-1], s₁⁻[k], p.λ)
            nfx₂[k] = numflux(s.u₂[k-1], s.u₂[k], p.f₂(s.u₂[k-1]), p.f₂(s.u₂[k]), s₂⁺[k-1], s₂⁻[k], p.λ)
            nfx₃[k] = numflux(s.u₃[k-1], s.u₃[k], p.f₃(s.u₃[k-1]), p.f₃(s.u₃[k]), s₃⁺[k-1], s₃⁻[k], p.λ)
        end

        # fluxes junction
        if p.cpl == TrafficFlowMaximization
            nfx₁[N] = v_R1
            nfx₂[N] = v_R2
            nfx₃[1] = v_L3
        else
            nfx₁[N] = numflux(s.u₁[end], u_R1, p.f₁(s.u₁[end]), v_R1, 0, 0, p.λ)
            nfx₂[N] = numflux(s.u₂[end], u_R2, p.f₂(s.u₂[end]), v_R2, 0, 0, p.λ)
            nfx₃[1] = numflux(u_L3, s.u₃[1], v_L3, p.f₃(s.u₃[1]), 0, 0, p.λ)
        end
        
        println("interface fluxes: street 1: ", nfx₁[[N-1, N]], ", street 2: ", nfx₂[[N-1, N]], ", street 3: ", nfx₃[[1, 2]])
        if abs(nfx₁[N] + nfx₂[N] - nfx₃[1]) > 1e-12 && ~warned
            println("Scheme not conservative")
            warned = true
        end
        
        # state update
        s.u₁ .-= dt/dx * diff(nfx₁)
        s.u₂ .-= dt/dx * diff(nfx₂)
        s.u₃ .-= dt/dx * diff(nfx₃)
        s.t += dt
        #readline()
    end
    s
end

    
function centraluncoupled(p::Problem11, s::NumSolution11, T::Real, sc::Scheme=FirstOrder, CFL::Real = 0.9)
    dx = (s.xc.stop - s.xc.start) / s.xc.lendiv
    if p.fixeddt > 0
        dt = p.fixeddt
    else
        dt = CFL * dx / p.λ
    end
    N = s.N

    nfx = zeros(2*N-1)

    s⁻ = zeros(2*N-2)
    s⁺ = zeros(2*N-2)
    sghost = (0, 0)
    
    while s.t<T

        ughost, vghost = ghostdata(s, p)
        uf = u(s)
        
        if sc == SecondOrder
            s⁻, s⁺ = charslopes(vcat(ughost[1], uf, ughost[2]), vcat(vghost[1], p.f₁.(uf), vghost[2]), p.λ)
            if p.bc == Periodic
                sghost = (s⁺[end], s⁻[1])
            end
        end

        if p.bc in [HomNeumann, Periodic]
            nfx[1] = numflux(ughost[1], s.u₁[1], vghost[1], p.f₁(s.u₁[1]), sghost[1],  s⁻[1], p.λ)
            nfx[end] = numflux(s.u₂[end], ughost[2], p.f₂(s.u₂[end]), vghost[2], s⁺[end], sghost[2], p.λ)
        end

        # fluxes
        for k=2:2*N-2
            nfx[k] = numflux(uf[k-1], uf[k], p.f₁(uf[k-1]), p.f₁(uf[k]), s⁺[k-1], s⁻[k], p.λ)
        end
                
        # state update
        #@show s⁻, s⁺
        #@show sum(diff(nfx))
        s.u₁ .-= dt/dx * diff(nfx[1:N])
        s.u₂ .-= dt/dx * diff(nfx[N:end])
        s.t += dt
        println("t=", round(s.t, digits=3), ", mass=", round(sum(abs.(u(s))) * dx, digits=3),
                ", TV=", round(sum(abs.(diff(u(s)))), digits=3))
    end
    s
end


end # module CentralNetworkScheme
