## Instructions for running the code in this folder

Point Julia's working directory to this folder using `cd("~\Linear_Hybrid_Latent_Dynamical_Models")`. Then:

1. Execute `include("Switch_Reactor_M1.jl")` to run filtering (using the bootstrap PF) on a linearised switching model of the CSTR. ONLY temperature is measured.
2. Execute `include("Switch_Reactor_M2.jl")` to run filtering (using the bootstrap PF) on a linearised switching model of the CSTR. Both temperature and concentration is measured.
3. Execute `include("Switch_Reactor_M1_RBPF.jl")` to run filtering (using the Rao-Blackwellised PF) on a linearised switching model of the CSTR. ONLY temperature is measured.
4. Execute `include("Switch_Reactor_M2_RBPF.jl")` to run filtering (using the Rao-Blackwellised PF) on a linearised switching model of the CSTR. Both temperature and concentration is measured.
