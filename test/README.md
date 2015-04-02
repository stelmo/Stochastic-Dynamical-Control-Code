## Running the tests

All testing is automated with Travis-CI but manually tests may also be performed by executing the following scripts:

1. `Reactor_test.jl` => Tests the integration routine (Runge-Kutta 4th order) of the coupled nonlinear CSTR system.
2. `HMM_test.jl` => Tests the inference algorithms of the Hidden Markov Model (HMM) section.
3. `LLDS_test.jl` => Tests the inference algorithms of the Linear Latent Dynamical Systems (LLDS) section.
4. `PF_test.jl` => Tests the SIR (bootstrap) Particle Filtering algorithm.
