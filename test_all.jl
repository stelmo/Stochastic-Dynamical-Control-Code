# Testing script

using Distributions

# Bring into scope all the required packages and functions.
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

include("./test/HMM_test.jl")
include("./test/LLDS_test.jl")
include("./test/PF_test.jl")
include("./test/Reactor_test.jl")
