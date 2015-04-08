# Bring into scope all the required packages and modules.
info("Loading all the required external modules.")
using Distributions
using PyPlot

info("Loading all the required internal modules.")
for func in ["Ellipse.jl",
            "HMM.jl",
            "Burglar.jl",
            "Reactor.jl",
            "LLDS.jl",
            "PF.jl",
            "SPF.jl",
            "RBPF.jl"]
  include(joinpath("src",func))
end
