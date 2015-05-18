# Stochastic Dynamical Control [![Build Status](https://travis-ci.org/stelmo/Stochastic-Dynamical-Control-Code.svg?branch=master)](https://travis-ci.org/stelmo/Stochastic-Dynamical-Control-Code)
Code for M. Eng. in Control Engineering

## Required Internal Modules:

**NB: Run the script `init_sys.jl` at startup to load all the required modules!**

## Required External Julia Packages:

1. Distributions
2. PyPlot
3. NLsolve
4. JuMP
5. Ipopt (including the optimisation library)
6. Mosek (including the optimisation library)
7. KernelDensity

## Running the Code:

Please see the folder `openloop_scenarions` or `closedloop_scenarions` to run simulations/demonstrations. The folders `openloop` and `closedloop` are meant as testing ground before scripts move into the scenario folders.

## Testing:

Please execute the script `test_all.jl` to perform the tests manually. This is done automatically with Travis-CI.
