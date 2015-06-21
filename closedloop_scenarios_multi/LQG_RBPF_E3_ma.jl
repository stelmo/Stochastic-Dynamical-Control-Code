# Control using multiple linear models and measuring both concentration and temperature

tend = 200
include("closedloop_params.jl") # load all the parameters and modules

# Get the three linear models about the nominal operating points
linsystems = Reactor.getNominalLinearSystems(h, cstr_model)

# Setup the RBPF
models, A = RBPF.setup_RBPF(linsystems, C2, Q, R2)
A = [0.99 0.01 0.00;
     0.01 0.98 0.01;
     0.00 0.01 0.99]
numModels = length(models) # number of linear models (will be 3)
nP = 500 # number of particles

init_state = linsystems[2].op

sguess =  RBPF.getInitialSwitches(init_state, linsystems) # prior switch distribution
particles = RBPF.init_RBPF(Categorical(sguess), init_state, init_state_covar, 2, nP)

maxtrack = zeros(length(linsystems), N) # keep track of the most likely model
switchtrack = zeros(length(linsystems), N) # keep track of the model/switch distribution

# Setup the controllers
setpoint = 0.998
H = [1.0 0.0]

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
horizon = 150
Anew, Bnew, bnew = Auxiliary.modelAverage(switchtrack[:, 1], models)
x_off, u_off = LQR.offset(Anew,Bnew, C2, H, ysp)
us[1] = MPC.mpc_lqr(rbpfmeans[:, 1] - bnew, horizon, Anew, Bnew, bnew, QQ, RR, x_off[1], u_off[1])
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
    Anew, Bnew, bnew = Auxiliary.modelAverage(switchtrack[:, t], models)
    x_off, u_off = LQR.offset(Anew,Bnew, C2, H, ysp)
    us[t] = MPC.mpc_lqr(rbpfmeans[:, t] - b, horizon, A, B, b, QQ, RR, x_off[1], u_off[1])
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
