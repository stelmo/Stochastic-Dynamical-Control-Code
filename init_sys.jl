# Bring into scope all the required packages and modules.
info("Loading all the required external modules.")
using PyPlot
using NLsolve
using Distributions

# You MUST run this from the root directory
# if nprocs() == 1
#   info("Adding additional processes.")
#   addprocs(CPU_CORES)
#   # addprocs(1)
# end
#
# @everywhere using Distributions
#
# info("Loading all the required internal modules.")
# @everywhere include(joinpath("src", "Ellipse.jl"))
# @everywhere include(joinpath("src", "HMM.jl"))
# @everywhere include(joinpath("src", "Burglar.jl"))
# @everywhere include(joinpath("src", "Reactor.jl"))
# @everywhere include(joinpath("src", "LLDS.jl"))
# @everywhere include(joinpath("src", "PF.jl"))
# @everywhere include(joinpath("src", "SPF.jl"))
# @everywhere include(joinpath("src", "RBPF.jl"))
# @everywhere include(joinpath("src", "LQR.jl"))
# @everywhere include(joinpath("src", "PSO.jl"))

for func in ["Ellipse.jl",
            "HMM.jl",
            "Burglar.jl",
            "Reactor.jl",
            "LLDS.jl",
            "PF.jl",
            "SPF.jl",
            "RBPF.jl",
            "LQR.jl",
            "PSO.jl"]
  include(joinpath("src",func))
end
