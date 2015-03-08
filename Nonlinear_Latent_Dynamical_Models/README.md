## Instructions for running the code in this folder

Point Julia's working directory to this folder using `cd("~\Nonlinear_Latent_Dynamical_Models")`. Then:

1. Execute `include("Nonlinear_Reactor.jl")` to run filtering and prediction on the CSTR reactor model. The nonlinear model is used to generate observations and a bootstrap particle filter is used to perform inference. Edit this script if you would like to change parameters.
