# Inference using one linear model measuring only temperature

tend = 50
include("openloop_params.jl") # load all the parameters and modules

init_state = [0.5, 400]

# Specify the linear model
linsystems = Reactor.getNominalLinearSystems(h, cstr_model)
opoint = 2 # which nominal model to use
A = linsystems[opoint].A
B = linsystems[opoint].B
b = linsystems[opoint].b

lin_cstr = LLDS.llds(A, B, C1, Q, R1) # KF object

# Plant initialisation
xs[:,1] = init_state
linxs[:, 1] = init_state - b

# Simulate plant
state_noise_dist = MvNormal(Q)
meas_noise_dist = Normal(0.0, sqrt(R1[1]))
ys1[1] = C1*xs[:, 1] + rand(meas_noise_dist) # measure from actual plant

# Filter setup
kfmeans = zeros(2, N)
kfcovars = zeros(2,2, N)
init_mean = init_state - b

# First time step
kfmeans[:, 1], kfcovars[:,:, 1] = LLDS.init_filter(init_mean, init_state_covar, ys1[1]-b[2], lin_cstr)

for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist)# actual plant
  ys1[t] = C1*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  linxs[:, t], temp = LLDS.step(linxs[:, t-1], us[t-1], lin_cstr)
  kfmeans[:, t], kfcovars[:,:, t] = LLDS.step_filter(kfmeans[:, t-1], kfcovars[:,:, t-1], us[t-1], ys1[t] - b[2], lin_cstr)
end

linxs = linxs .+ b
kfmeans = kfmeans .+ b

# Plot results

Results.plotEllipses(ts, xs, kfmeans, kfcovars, "Kalman Filter", "upper right")

Results.plotTracking(ts, xs, ys1, kfmeans, us, 1)

avediff = Results.calcError(xs, kfmeans)
