# Closed Loop Scenarios

This folder contains stand alone scripts for illustration purposes. They are derived from the scripts in the folder `closedloop`.

**It is necessary to run the script `init_sys.jl` in the root directory to load all the external and internal libraries.**

All the estimated run times are relative to my laptop. The technical specifications of it are: 4 GB RAM, 2.20GHz Intel Core i7 CPU (26270QM), Windows 7.

## Single Model LQG

## Switching Model LQG

## Linear system model controlled with a linear MPC

1. `lin_mod_kf_lin_mpc_mean.jl` => Linear plant model controlled with a linear MPC. Plant and measurement noise is Gaussian with known distributions. A Kalman Filter is used for state estimation. The control objective is to keep the system at the unsteady operating point using the minimum controller effort (standard QP MPC objective function). Standard deterministic constraints. Takes about 9 seconds to simulate 20 minutes.

2. `lin_mod_kf_lin_mpc_var_conf_x.jl` => exactly the same system as before but with the additional stochastic constraint. The suffix `_x.jl` indicates the confidence of the constraint based on the Gaussian assumption i.e. it is read off from a Chi Squared Distribution table (with 2 DOF in this case). Takes about 10 seconds to simulate 20 minutes.

3. `lin_mod_pf_lin_mpc_mean.jl` => The same as (1) but using a Particle Filter for state estimation. Using 2000 particles it takes about 34 seconds to simulate 20 minutes.

4. `lin_mod_pf_lin_mpc_var_conf_99_graphical.jl` => Similar to (2) but we additionally test the assumption that the posterior distributions are Gaussian. The density function of the posterior is estimated from the PF samples and the contours are plotted over time intervals. This allows one to qualitatively inspect the Gaussian nature of the pdfs. Takes about 9 minutes using 70000 particles and simulating 50 min.

5. `lin_mod_pf_lin_mpc_var_conf_99_numeric.jl` => Similar to (2) but we additionally test the assumption that the posterior distributions are Gaussian. The Kullback-Leibler Divergence of the sampled posterior is compared to the Gaussian simplification thereof. The baseline for this test can be found in `openloop_scenarios/KL_PF_lin_mod_baseline.jl`. Note that this test takes approximately 19 min to complete if the total simulated time is 20 min and the number of particles is 20000.

## Nonlinear system model controlled with a linear MPC

1. `nonlin_mod_kf_lin_mpc_mean.jl` => Nonlinear plant model controlled with a linear MPC. Plant and measurement noise is Gaussian with known distributions. A Kalman Filter is used for state estimation. The control objective is to keep the system at the unsteady operating point using the minimum controller input (standard QP MPC objective function). Standard deterministic constraints are used.

2. `nonlin_mod_kf_lin_mpc_var_conf_x.jl` => Exactly the same system as before but with the additional stochastic constraint. The suffix `_x.jl` indicates the confidence of the constraint based on the Gaussian assumption i.e. it is read off from a Chi Squared Distribution table (with 2 DOF in this case). Takes about 24 seconds to simulate 50 minutes.

3. `lin_mod_pf_lin_mpc_mean.jl` => The same as (1) but using a Particle Filter for state estimation. Using 7000 particles it takes about 63 seconds to simulate 50 minutes.

4. `nonlin_mod_pf_lin_mpc_var_conf_x.jl` => Exactly the same system as (3) but with the additional stochastic constraint. The suffix `_x.jl` indicates the confidence of the constraint based on the Gaussian assumption. Takes about 70 seconds to simulate 50 minutes using 7000 particles.

5. `lin_mod_pf_lin_mpc_var_conf_99_graphical.jl` => Similar to (2) but we additionally test the assumption that the posterior distributions are Gaussian. The density function of the posterior is estimated from the PF samples and the contours are plotted over time intervals. This allows one to qualitatively inspect the Gaussian nature of the pdfs. Takes about 10 minutes using 70000 particles and simulating 50 min.

6. `lin_mod_pf_lin_mpc_var_conf_99_numeric.jl` => Similar to (2) but we additionally test the assumption that the posterior distributions are Gaussian. The Kullback-Leibler Divergence of the sampled posterior is compared to the Gaussian simplification thereof. The baseline for this test can be found in `openloop_scenarios/KL_PF_lin_mod_baseline.jl`. Note that this test takes approximately 19 min to complete if the total simulated time is 20 min and the number of particles is 20000.


## Switching Model Linear MPC
