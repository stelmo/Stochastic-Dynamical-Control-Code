# Test the Particle Filter.
# We conduct the tests by comparing the posterior
# filtered densities to the analytic Kalman Filter
# solution as calculated by the functions in the
# Linear_Latent_Dynamical_Models folder.

import PF
reload("PF.jl")

nP = 10 #number of particles.
init_state_mean = zeros(2) # initial state mean
init_state_covar = eye(2) # initial state covariance
dist = MvNormal(init_state_mean, init_state_covar) # prior distribution
particles = PF.init_PF(dist, nP, 1) # initialise the particles
