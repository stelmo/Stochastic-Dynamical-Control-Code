# Stochastic Dynamical Control [![Build Status](https://travis-ci.org/stelmo/Stochastic-Dynamical-Control-Code.svg?branch=master)](https://travis-ci.org/stelmo/Stochastic-Dynamical-Control-Code)
Code for M. Eng. in Control Engineering

## Required Julia Packages:

1. Distributions
2. PyPlot
3. Gadfly

## Running the Code

Please see the README inside each folder for details about running the code in each section.

## Testing
All testing is automated with Travis-CI but manually tests may also be performed by executing the following scripts:

1. Inside the `Hidden_Markol_Models` folder: execute the script `HMM_tests.jl`

2. Inside the `CSTR_Model` folder: execute the script `Reactor_test.jl`

3. Inside the `Linear_Latent_Dynamic_Models` folder: execute the script `LLDS_tests.jl`

4. Inside the `Noninear_Latent_Dynamical_Models` folder: execute the script `Nonlinear_test.jl`

5. Inside the `Linear_Hybrid_Latent_Dynamical_Models` folder: execute the script `Linear_Switch_test.jl`

6. Inside the `Nonlinear_Hybrid_Latent_Dynamical_Models` folder: execute the script `Nonlinear_Switch_test.jl`
