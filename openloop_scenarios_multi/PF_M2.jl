# PF inference using the full nonlinear model

tend = 300
include("openloop_params.jl") # load all the parameters and modules
srand(852)
init_state = [0.55; 450] # initial state

f(x, u, w) = Reactor.run_reactor(x, u, h, cstr_model) + w
g(x) = C2*x # state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 200 # number of particles.
prior_dist = MvNormal(init_state, init_state_covar) # prior distribution
particles = PF.init_PF(prior_dist, nP, 2) # initialise the particles
state_noise_dist = MvNormal(Q) # state distribution
meas_noise_dist = MvNormal(R2) # measurement distribution

# Time step 1
xs[:, 1] = init_state
xsnofix[:, 1] = init_state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys2[:, 1], meas_noise_dist, cstr_pf)
pfmeans[:,1], pfcovars[:,:,1] = PF.getStats(particles)

# Loop through the rest of time
for t=2:N
  random_element = rand(state_noise_dist)
  if ts[t] < 50.0
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
  else
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model_broken) + random_element
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
  end
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  PF.filter!(particles, us[t-1], ys2[:, t], state_noise_dist, meas_noise_dist, cstr_pf)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles)
end


# Plot results
Results.plotTrackingBreak(ts, xs, xsnofix, ys2, pfmeans, 2)
Results.calcError(xs, pfmeans)
