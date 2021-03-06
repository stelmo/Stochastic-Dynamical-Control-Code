# Control using multiple linear models and measuring both concentration and temperature

tend = 200
include("closedloop_params.jl") # load all the parameters and modules

# Divide state space into sectors: n by m
nX = 2 # rows
nY = 2 # cols
xspace = [0.0, 1.0]
yspace = [250, 550]

linsystems = Reactor.getLinearSystems(nX, nY, xspace, yspace, h, cstr_model)

models, A = RBPF.setup_RBPF(linsystems, C2, Q, R2)
A = [0.98 0.00 0.00 0.00 0.00 0.01 0.00;
     0.00 0.98 0.00 0.01 0.01 0.01 0.00;
     0.01 0.00 0.98 0.00 0.00 0.01 0.01;
     0.00 0.00 0.00 0.98 0.00 0.01 0.00;
     0.00 0.01 0.00 0.00 0.99 0.00 0.00;
     0.01 0.01 0.01 0.01 0.00 0.96 0.00;
     0.00 0.00 0.01 0.00 0.00 0.00 0.99]
numModels = length(models)

nP = 500 # number of particles
init_state = linsystems[6].op

sguess =  RBPF.getInitialSwitches(init_state, linsystems) # prior switch distribution
particles = RBPF.init_RBPF(Categorical(sguess), init_state, init_state_covar, 2, nP)

maxtrack = zeros(length(linsystems), N)
switchtrack = zeros(length(linsystems), N)

# # Setup the controllers
setpoint = 0.9 #linsystems[3].op[1]
H = [1.0 0.0]
controllers = Array(LQR.controller, length(models))
for k=1:length(models)
  ysp = setpoint - models[k].b[1] # set point is set here
  x_off, u_off = LQR.offset(models[k].A,models[k].B, C2, H, ysp)
  K = LQR.lqr(models[k].A, models[k].B, QQ, RR)
  controllers[k] = LQR.controller(K, x_off, u_off)
end

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

# Controller Input
ind = indmax(maxtrack[:, 1]) # use this model and controller
horizon = 150
us[1] = MPC.mpc_lqr(rbpfmeans[:, 1] - models[ind].b, horizon, models[ind].A, models[ind].B, models[ind].b, QQ, RR, controllers[ind].x_off[1], controllers[ind].u_off[1])

#Loop through the rest of time
tic()
for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  RBPF.filter!(particles, us[t-1], ys2[:, t], models, A)
  rbpfmeans[:,t], rbpfcovars[:,:, t] = RBPF.getAveStats(particles)

  for k=1:numModels
    switchtrack[k, t] = sum(particles.ws[find((x)->x==k, particles.ss)])
  end
  maxtrack[:, t] = RBPF.getMaxTrack(particles, numModels)

  # Controller Input
  if t%10 == 0
    ind = indmax(maxtrack[:, t]) # use this model and controller
    us[t] = MPC.mpc_lqr(rbpfmeans[:, t] - models[ind].b, horizon, models[ind].A, models[ind].B, models[ind].b, QQ, RR, controllers[ind].x_off[1], controllers[ind].u_off[1])
  else
    us[t] = us[t-1]
  end
end
toc()

# Plot results
Results.plotStateSpaceSwitch(linsystems, xs)
Results.plotSwitchSelection(numModels, maxtrack, ts, false)
Results.plotSwitchSelection(numModels, switchtrack, ts, true)
Results.plotTracking(ts, xs, ys2, rbpfmeans, us, 2, setpoint)
