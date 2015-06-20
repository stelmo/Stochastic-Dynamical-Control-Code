# Implement the augmented switching dynamical system

tend = 150
include("openloop_params.jl") # load all the parameters and modules

init_state = [0.5; 400] # initial state

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

nP = 500
sguess =  RBPF.getInitialSwitches(init_state, linsystems)
particles = RBPF.init_RBPF(Categorical(sguess), init_state, init_state_covar, 2, nP)

switchtrack = zeros(length(linsystems), N)
maxtrack = zeros(length(linsystems), N)
smoothedtrack = zeros(length(linsystems), N)

state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(R2)

xs[:,1] = init_state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant

RBPF.init_filter!(particles, 0.0, ys2[:, 1], models)
rbpfmeans[:,1], rbpfcovars[:,:, 1] = RBPF.getAveStats(particles)


for k=1:length(linsystems)
  switchtrack[k, 1] = sum(particles.ws[find((x)->x==k, particles.ss)])
end
maxtrack[:, 1] = RBPF.getMaxTrack(particles, numModels)
smoothedtrack[:, 1] = RBPF.smoothedTrack(numModels, switchtrack, 1, 10)

# Loop through the rest of time
tic()
for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  RBPF.filter!(particles, us[t-1], ys2[:, t], models, A)
  rbpfmeans[:, t], rbpfcovars[:,:, t] = RBPF.getAveStats(particles)

  for k=1:length(linsystems)
    switchtrack[k, t] = sum(particles.ws[find((x)->x==k, particles.ss)])
  end
  maxtrack[:, t] = RBPF.getMaxTrack(particles, numModels)
  smoothedtrack[:, t] = RBPF.smoothedTrack(numModels, switchtrack, t, 10)

end
toc()

# Plot results
Results.plotStateSpaceSwitch(linsystems, xs)
Results.plotSwitchSelection(numModels, switchtrack, ts, true)
# Results.plotSwitchSelection(numModels, maxtrack, ts, false)
# Results.plotSwitchSelection(numModels, smoothedtrack, ts, false)
Results.plotTracking(ts, xs, ys2, rbpfmeans, us, 2)
