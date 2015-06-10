# Open Loop Scenarios

This folder contains stand alone scripts for illustration purposes. They are derived from the scripts in the folder `closedloop`.

**It is necessary to run the script `init_sys.jl` in the root directory to load all the external and internal libraries.**

All the estimated run times are relative to my laptop. The technical specifications of it are: 4 GB RAM, 2.20GHz Intel Core i7 CPU (26270QM), Windows 7.

1. `KL_PF_lin_mod_baseline_graphical.jl` => Plots the pdf of the posterior distribution as a function of time using the PF. Since the model is linear and the noise Gaussian this should be approximately Gaussian. Zero inputs are used. Using 70000 particles and simulating 20 minutes takes about 3 minutes.
2. `KL_PF_lin_mod_baseline_numerical.jl` => The same as (1) except that the Kullback Leibler divergence of the estimated pdf and the Gausssian approximation is plotted. Zero inputs are used. Using 70000 particles and simulating 20 minutes this takes about 18 minutes.

### Burglar Example: Analytic Hidden Markov Model Inference

1. Execute `HMM_burglar_example.jl` to generate the HMM figures.  

### Linear Reactor: Kalman Filtering

In this section the uncertain, noisy nonlinear CSTR model generates the observations and the unstable linear model is used for inference.

1. Execute `Linear_Reactor_M1.jl` to perform Kalman Filtering by measuring only temperature.
2. Execute `Linear_Reactor_M2.jl` to perform Kalman Filtering by measuring both temperature and concentration.

### Nonlinear Reactor: Particle Filtering

In this section the uncertain, noisy nonlinear CSTR model generates the observations and the nonlinear model is used to perform inference using Particle techniques.

1. Execute `Nonlinear_Reactor_M1.jl` to perform Particle Filtering by measuring only temperature.
2. Execute `Nonlinear_Reactor_M2.jl` to perform Particle Filtering by measuring both temperature and concentration.
3. Execute `Nonlinear_Linear_Reactor.jl` to compare the Particle and Kalman Filter using both measurements and the unstable linear model for inference.
