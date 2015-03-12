## Instructions for running the code in this folder

Point Julia's working directory to this folder using `cd("~\Nonlinear_Latent_Dynamical_Models")`. Then:

1. Execute `include("PF_test.jl")` to run the test.
2. Execute `include("Nonlinear_Reactor_M1.jl")` to perform filtering on the CSTR reactor model. The nonlinear model is used to generate observations and a bootstrap particle filter is used to perform inference using the nonlinear model. We ONLY measure temperature here. Edit this script if you would like to change parameters.
3. Execute `include("Nonlinear_Reactor_M2.jl")` to perform filtering on the CSTR reactor model. The nonlinear model is used to generate observations and a bootstrap particle filter is used to perform inference using the nonlinear model. We measure BOTH temperature and concentration here. Edit this script if you would like to change parameters.
4. Execute `include("PF_KF_Comparison.jl")` to perform filtering on the CSTR reactor model. The nonlinear model is used to generate observations and a bootstrap particle filter is compared to a Kalman Filter using a linear model of the process. Edit this script if you would like to change parameters.
