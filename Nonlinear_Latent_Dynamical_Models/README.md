## Instructions for running the code in this folder

Point Julia's working directory to this folder using `cd("~\Nonlinear_Latent_Dynamical_Models")`. Then:

1. Execute `include("Nonlinear_Reactor.jl")` to perform filtering on the CSTR reactor model. The nonlinear model is used to generate observations and a bootstrap particle filter is used to perform inference using the nonlinear model. Edit this script if you would like to change parameters.
2. Execute `include("Linear_Model.jl")` to perform filtering on the CSTR reactor model. The nonlinear model is used to generate observations and a bootstrap particle filter is used to perform inference using the linear model. Edit this script if you would like to change parameters.
