# Open Loop Scenarios

This folder contains stand alone scripts for illustration purposes. They are derived from the scripts in the folder `closedloop`.

**It is necessary to run the script `init_sys.jl` in the root directory to load all the external and internal libraries.**

All the estimated run times are relative to my laptop. The technical specifications of it are: 4 GB RAM, 2.20GHz Intel Core i7 CPU (26270QM), Windows 7.

1. `KL_PF_lin_mod_baseline_graphical.jl` => Plots the pdf of the posterior distribution as a function of time using the PF. Since the model is linear and the noise Gaussian this should be approximately Gaussian. Zero inputs are used. Using 2000 particles this takes about [] minutes.
2. `KL_PF_lin_mod_baseline_numerical.jl` => The same as (1) except that the Kullback Leibler divergence of the estimated pdf and the Gausssian approximation is plotted. Zero inputs are used. Using 2000 particles this takes about [] minutes. 
