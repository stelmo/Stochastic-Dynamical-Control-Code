# Julia Scripts
This folder contains all the scripts used to generate the figures as shown in my masters.

### Burglar Example: Analytic Hidden Markov Model Inference

1. Execute `HMM_burglar_example.jl` to generate the HMM figures.  


### Linear Reactor: Kalman Filtering and Prediction

In this section the uncertain, noisy nonlinear CSTR model generates the observations and the unstable linear model is used for inference.

1. Execute `Linear_Reactor_M1.jl` to perform Kalman Filtering and Prediction by measuring only temperature.
2. Execute `Linear_Reactor_M2.jl` to perform Kalman Filtering and Prediction by measuring both temperature and concentration.

### Nonlinear Reactor: Particle Filtering and Prediction

In this section the uncertain, noisy nonlinear CSTR model generates the observations and the nonlinear model is used to perform inference using Particle techniques.

1. Execute `Nonlinear_Reactor_M1.jl` to perform Particle Filtering and Prediction by measuring only temperature.
2. Execute `Nonlinear_Reactor_M2.jl` to perform Particle Filtering and Prediction by measuring both temperature and concentration.
3. Execute `Nonlinear_Linear_Reactor.jl` to compare the Particle and Kalman Filter using both measurements and the unstable linear model for inference.

### Switching Linear Reactor: Rao-Blackwellised Particle Filter

In this section the uncertain, noisy nonlinear CSTR model generates the observations. The state space is divided into regions and the center of each region is used to create a linear model. These linear models are the constitutive models used in the switching framework. Rao-Blackwellisation is used to perform inference.

1. Execute `Switching_Linear_Reactor_M1.jl` to perform Filtering and Prediction by measuring only temperature.
2. Execute `Switching_Linear_Reactor_M2.jl` to perform Filtering and Prediction by measuring both temperature and concentration.

### Switching Nonlinear Reactor: Switching Particle Filter

In this section uncertain, noisy nonlinear CSTR models generate the observations. We assume that we have two nonlinear models describing the system. The system initially uses the first model but at a certain time "breaks". The underlying process is then best described by the second nonlinear model. These to nonlinear models are the constitutive models for the underlying switching framework.

1. Execute `Switching_Nonlinear_Reactor_M1.jl` to perform Filtering and Prediction by measuring only temperature.
2. Execute `Switching_Nonlinear_Reactor_M2.jl` to perform Filtering and Prediction by measuring both temperature and concentration.
