# Bring into scope all the required packages and modules.
info("Loading all the required external modules.")
using PyPlot
using NLsolve
using Distributions

for func in ["Ellipse.jl",
            "HMM.jl",
            "Burglar.jl",
            "Reactor.jl",
            "LLDS.jl",
            "PF.jl",
            "SPF.jl",
            "RBPF.jl",
            "LQR.jl",
            "PSO.jl",
            "Results.jl",
            "MPC.jl"]
  include(joinpath("src",func))
end
