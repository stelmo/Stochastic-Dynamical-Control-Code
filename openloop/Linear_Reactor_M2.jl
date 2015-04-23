# Inference using the linear reactor model measuring both concentration and temperature

include("../params.jl") # load all the parameters and modules

init_state = [0.50, 400]

# Specify the linear model
linsystems = Reactor.getNominalLinearSystems(h, cstr_model)
opoint = 2 # which nominal operating point to use

A = linsystems[opoint].A
B = linsystems[opoint].B
b = linsystems[opoint].b

lin_cstr = LLDS.llds(A, B, C2, Q, R2) # KF object

# Plant initialisation
xs[:,1] = init_state
linxs[:, 1] = init_state - b

# Simulate plant
state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(lin_cstr.R)
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measure from actual plant

# Filter initialisation
kfmeans = zeros(2, N)
kfcovars = zeros(2,2, N)
init_mean = init_state - b
init_covar = eye(2)
init_covar[1] = 1e-3
init_covar[4] = 4.

# First time step
kfmeans[:, 1], kfcovars[:,:, 1] = LLDS.init_filter(init_mean, init_covar, ys2[:, 1]-b, lin_cstr)

for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  linxs[:, t], temp = LLDS.step(linxs[:, t-1], us[t-1], lin_cstr)
  kfmeans[:, t], kfcovars[:,:, t] = LLDS.step_filter(kfmeans[:, t-1], kfcovars[:,:, t-1], us[t-1], ys2[:, t]-b, lin_cstr)
end

linxs = linxs .+ b
kfmeans = kfmeans .+ b

# Plot results
Results.plotEllipses(ts, xs, ys2, kfmeans, kfcovars)

Results.plotTracking(ts, xs, ys2, kfmeans, us, 2)
