# Scenarios

This folder contains stand alone scripts for illustration purposes. They are derived from the scripts in the folder `closedloop`.

**It is necessary to run the script `init_sys.jl` in the root directory to load all the external and internal libraries.**

## Single Model LQG

## Switching Model LQG



## Single Model Linear MPC

1. `lin_mod_lin_mpc_gauss_mean.jl` => Linear plant model controlled with a linear MPC. Plant and measurement noise is Gaussian with known distributions. The control objective is to keep the system at the unsteady operating point using the minimum controller effort (standard QP MPC objective function). Standard deterministic constraints.

2. `lin_mod_lin_mpc_gauss_var_conf_x.jl` => exactly the same system as before but with the additional stochastic constraint. The suffix `_x.jl` indicates the confidence of the constraint based on the Gaussian assumption i.e. it is read off from a Chi Squared Distribution table (with 2 DOF in this case).





## Switching Model Linear MPC
