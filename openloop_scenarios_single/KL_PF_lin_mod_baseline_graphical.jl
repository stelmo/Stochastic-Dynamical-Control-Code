# PF inference using the full nonlinear model. Illustration of Kullback Leibler divergence.
# It should be close to zero here.

tend = 20
include("openloop_params.jl") # load all the parameters and modules

# Get the linear model
linsystems = Reactor.getNominalLinearSystems(h, cstr_model) # cstr_model comes from params.jl
opoint = 2 # the specific operating point we are going to use for control

init_state = [0.5, 400] # random initial point near operating point

# Set the state space model
A = linsystems[opoint].A
B = linsystems[opoint].B
b = linsystems[opoint].b # offset from the origin

f(x, u, w) = A*x + B*u + w
g(x) = C2*x # state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 70000 # number of particles.
prior_dist = MvNormal(init_state-b, init_state_covar) # prior distribution
particles = PF.init_PF(prior_dist, nP, 2) # initialise the particles
state_noise_dist = MvNormal(Q) # state distribution
meas_noise_dist = MvNormal(R2) # measurement distribution

# Time step 1
xs[:,1] = init_state-b
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys2[:, 1], meas_noise_dist, cstr_pf)
pfmeans[:,1], pfcovars[:,:,1] = PF.getStats(particles)

temp_states = zeros(2, nP)
# Loop through the rest of time
tic()
for t=2:N
  xs[:, t] = A*xs[:, t-1] + B*us[t-1] + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  PF.filter!(particles, us[t-1], ys2[:, t], state_noise_dist, meas_noise_dist, cstr_pf)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles)

  if ts[t] in [0.0:3.0:tend]
    Auxiliary.showEstimatedDensity(particles.x, particles.w, temp_states)
  end
end
toc()
xs = xs .+ b
pfmeans = pfmeans .+ b
ys2 = ys2 .+ b

# Plot results
Results.plotTracking(ts, xs, ys2, pfmeans, us, 2)
