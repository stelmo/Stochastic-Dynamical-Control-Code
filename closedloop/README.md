# Closed Loop Systems

This folder contains all the scripts used to generate the closed loop figures as shown in my masters.

### Linear Quadratic Gaussian Control

1. `LQG_KF.jl` => LQG control using the Kalman Filter as a state estimator.
2. `LQG_RBPF.jl` => LQG control switching controller models based on the Rao-Blackwellised Particle Filter state estimate.
3. `LQG_SPF.jl` => LQG control switching controller models based on the Switching Particle filter. In this scenario the plant breaks.
4. `LQG_SPF_COMP.jl` => LQG control showing the benefit of using the SPF vs. not switching models when the plant breaks.

### Linear MPC

1. `MPC_KF_MEAN.jl` => MPC using the KF as the state estimator with traditional deterministic constraints.
2. `MPC_KF_VAR.jl` => MPC using the KF as the state estimator with stochastic (chance) constraints.
3. `MPC_RBPF_MEAN.jl` => Switching MPC using the Rao Blackwellised Particle Filter to select which model to use. Deterministic constraints.
4. `MPC_RBPF_VAR.jl` => Switching MPC using the Rao Blackwellised Particle Filter to select which model to use. Stochastic constraints.
5. `MPC_SPF_MEAN.jl` => Switching MPC using the Switching Particle Filter to select which model to use. Deterministic constraints. In this scenario the plant breaks.
6. `MPC_SPF_VAR.jl` => Switching MPC using the Switching Particle Filter to select which model to use. Stochastic constraints. In this scenario the plant breaks.

### Nonlinear MPC

1.
2.
