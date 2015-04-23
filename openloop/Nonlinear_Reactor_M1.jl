# PF Inference using the full nonlinear model

include("../params.jl") # load all the parameters and modules

init_state = [0.5; 400] # initial state

f(x, u, w) = Reactor.run_reactor(x, u, h, cstr_model) + w # transition function
g(x) = C1*x # state observation

cstr_pf = PF.Model(f,g) # PF object

# Initialise the PF
nP = 100 # number of particles.
prior_dist = MvNormal(init_state, init_state_covar) # prior distribution
particles = PF.init_PF(prior_dist, nP, 2) # initialise the particles

state_noise_dist = MvNormal(Q) # state distribution
meas_noise_dist = MvNormal(eye(1)*R1) # measurement distribution

# Time step 1
xs[:, 1] = init_state
ys1[1] = [0.0 1.0]*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys1[1], meas_noise_dist, cstr_pf)
pfmeans[:,1], pfcovars[:,:,1] = PF.getStats(particles)
# Loop through the rest of time
for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys1[t] = C1*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  PF.filter!(particles, us[t-1], ys1[t], state_noise_dist, meas_noise_dist, cstr_pf)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles)
end

# Plot results
Results.plotEllipses(ts, xs, ys1, pfmeans, pfcovars)

Results.plotTracking(ts, xs, ys1, pfmeans, us, 1)
