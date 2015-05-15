# Scenarios

This folder contains stand alone scripts for illustration purposes. They are derived from the scripts in the folder `closedloop`.

**It is necessary to run the script `init_sys.jl` in the root directory to load all the external and internal libraries.**

## Single Model LQG

## Switching Model LQG



## Single Model Linear MPC

1. `lin_mod_lin_mpc_gauss.jl` => Linear plant model controlled with a linear MPC. Plant and measurement noise is Gaussian with known distributions. The control objective is to keep the system at the unsteady operating point using the minimum controller effort (standard QP MPC objective function).

## Switching Model Linear MPC
