# Control using the nonlinear reactor model measuring both concentration and temperature
# NOTE: this is very slow: 1 min simulated costs about 3 min in real life (parameter dependent)
include("../params.jl") # load all the parameters and modules

# Set up the PF
f(x, u, w) = Reactor.run_reactor(x, u, h, cstr_model) + w # transistion function
g(x) = C2*x # measurement function
pf_cstr = PF.Model(f,g) # PF object

# Simulation setup
init_state = [0.5; 400] # initial state

# Initialise the PF
nP = 20 # number of particles.
prior_dist = MvNormal(init_state, init_state_covar) # prior distribution
particles = PF.init_PF(prior_dist, nP, 2) # initialise the particles

state_noise_dist = MvNormal(Q) # state distribution
meas_noise_dist = MvNormal(R2) # measurement distribution

# Controller setup
y_ca = 0.48 # concentration set point (0.48934869384882873)
offres = PSO.offset(y_ca, cstr_model)
ys2p = [y_ca, offres.zero[1]]
usp = offres.zero[2]

# PSO setup
skip = 7 # hold the controller output for this long
predictionHorizon = 10 # number of controller moves predicted/optimised (control horizon)
swarmSize = 50 # number of particles in the swarm
swarmOptRepeat = 100 # high is better but it slows everything down

# Time step 1
xs[:,1] = init_state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys2[:, 1], meas_noise_dist, pf_cstr)
pfmeans[:,1], pfcovars[:,:,1] = PF.getStats(particles)
# Controlle action
swarm, sol = PSO.initswarm(swarmSize, predictionHorizon, -1000.0, 1000.0, particles, ys2p, usp, QQ, RR, state_noise_dist, cstr_model, skip, h)
us[1] = PSO.optimise!(swarm, sol, particles, ys2p, usp, QQ, RR, state_noise_dist, cstr_model, skip, h, swarmOptRepeat)

for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  PF.filter!(particles, us[t-1], ys2[:, t], state_noise_dist, meas_noise_dist, pf_cstr)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles)

  # Controller action
  if t%skip == 0
    swarm, sol = PSO.initswarm(swarmSize, predictionHorizon, -1000.0, 1000.0, particles, ys2p, usp, QQ, RR, state_noise_dist, cstr_model, skip, h)
    us[t] = PSO.optimise!(swarm, sol, particles, ys2p, usp, QQ, RR, state_noise_dist, cstr_model, skip, h, swarmOptRepeat)
    println("Time: ", ts[t])
  end
  if t%skip != 0
    us[t] = us[t-1]
  end
  # println("Controller Input: ", us[t-1])

end

#  Plot Results
Results.plotTracking(ts, xs, ys2, pfmeans, us, 2)
