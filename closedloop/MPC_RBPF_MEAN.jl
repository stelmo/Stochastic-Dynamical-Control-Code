# Control using multiple linear models and measuring both concentration and temperature

include("../params.jl") # load all the parameters and modules

# Get the three linear models about the nominal operating points
linsystems = Reactor.getNominalLinearSystems(h, cstr_model) # cstr_model comes from params.jl
opoint = 2 # the specific operating point we are going to use for control

init_state = [0.5, 450] # random initial point near operating point

# Set point
ysp = linsystems[1].op[1] # Low concentration
# ysp = linsystems[2].op[1] # Medium concentration

# Setup the RBPF
models, A = RBPF.setup_RBPF(linsystems, C2, Q, R2)
numModels = length(models) # number of linear models (will be 3)
nP = 1500 # number of particles

sguess =  RBPF.getInitialSwitches(init_state, linsystems) # prior switch distribution
particles = RBPF.init_RBPF(Categorical(sguess), init_state, init_state_covar, 2, nP)

maxtrack = zeros(length(linsystems), N) # keep track of the most likely model
switchtrack = zeros(length(linsystems), N) # keep track of the model/switch distribution
smoothedtrack = zeros(length(linsystems), N)

# Setup MPC
horizon = 150
# add state constraints
# 540 @ 0.0 and 440 @ 0.8
aline = 125.0 # slope of constraint line ax + by + c = 0
cline = -540.0 # negative of the y axis intercept
bline = 1.0

state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(R2)
xs[:,1] = init_state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant

RBPF.init_filter!(particles, 0.0, ys2[:, 1], models)
rbpfmeans[:,1], rbpfcovars[:,:, 1] = RBPF.getAveStats(particles)

for k=1:numModels
  switchtrack[k, 1] = sum(particles.ws[find((x)->x==k, particles.ss)])
end
maxtrack[:, 1] = RBPF.getMaxTrack(particles, numModels)
smoothedtrack[:, 1] = RBPF.smoothedTrack(numModels, switchtrack, 1, 10)

# Controller Input
ind = indmax(smoothedtrack[:, 1]) # use this model and controller
yspfix = ysp - linsystems[ind].b[1]
us[1] = MPC.mpc_mean(rbpfmeans[:, 1]-linsystems[ind].b[1], horizon, linsystems[ind].A, linsystems[ind].B, linsystems[ind].b, aline, bline, cline, QQ, RR, yspfix, 20000.0, true) # get the controller input

#Loop through the rest of time
for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  RBPF.filter!(particles, us[t-1], ys2[:, t], models, A)
  rbpfmeans[:,t], rbpfcovars[:,:, t] = RBPF.getAveStats(particles)

  for k=1:numModels
    switchtrack[k, t] = sum(particles.ws[find((x)->x==k, particles.ss)])
  end
  maxtrack[:, t] = RBPF.getMaxTrack(particles, numModels)
  smoothedtrack[:, t] = RBPF.smoothedTrack(numModels, switchtrack, t, 20)

  # Controller Input
  ind = indmax(smoothedtrack[:, t]) # use this model and controller
  yspfix = ysp - linsystems[ind].b[1]
  us[t] = MPC.mpc_mean(rbpfmeans[:, t]-linsystems[ind].b, horizon, linsystems[ind].A, linsystems[ind].B, linsystems[ind].b, aline, bline, cline, QQ, RR, yspfix, 20000.0, true)# get the controller input
end

# Plot results
Results.plotStateSpaceSwitch(linsystems, xs)

Results.plotSwitchSelection(numModels, maxtrack, ts, false)

Results.plotSwitchSelection(numModels, smoothedtrack, ts, false)

Results.plotTracking(ts, xs, ys2, rbpfmeans, us, 2, ysp)

Results.plotEllipses(ts, xs, rbpfmeans, rbpfcovars, "MPC", [aline, cline], linsystems[1].op, true)
