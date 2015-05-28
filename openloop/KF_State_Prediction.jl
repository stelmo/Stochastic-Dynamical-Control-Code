# Controller using the linear reactor model measuring both concentration and temperature.

include("../params.jl") # load all the parameters and modules

# Get the linear model
linsystems = Reactor.getNominalLinearSystems(h, cstr_model) # cstr_model comes from params.jl
opoint = 2 # the specific operating point we are going to use for control

init_state = [0.5, 450] # random initial point near operating point

# Set the state space model
A = linsystems[opoint].A
B = linsystems[opoint].B
b = linsystems[opoint].b # offset from the origin

# Set up the KF
kf_cstr = LLDS.llds(A, B, C2, Q, R2) # set up the KF object (measuring both states)
state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(R2)

# Setup the prediction
pskip = 1
pnum = 150
pstart = 5
pmeans = zeros(2, pnum)
pcovars = zeros(2,2, pnum)

# First time step of the simulation
xs[:,1] = init_state # set simulation starting point to the random initial state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measure from actual plant
kfmeans[:, 1], kfcovars[:,:, 1] = LLDS.init_filter(init_state-b, init_state_covar, ys2[:, 1]-b, kf_cstr) # filter

for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measure from actual plant
  kfmeans[:, t], kfcovars[:,:, t] = LLDS.step_filter(kfmeans[:, t-1], kfcovars[:,:, t-1], us[t-1], ys2[:, t]-b, kf_cstr)

  # prediction
  if t == pstart
    pmeans, pcovars = LLDS.predict_hidden(kfmeans[:, t], kfcovars[:,:, t], us[t:1:t+pnum-1], kf_cstr)
  end
end
kfmeans = kfmeans .+ b
pmeans = pmeans .+ b

# Plot the results
# Results.plotTracking(ts, xs, ys2, kfmeans, us, 2)

# Results.plotEllipses(ts, xs, kfmeans, kfcovars, "Kalman Filter")

Results.plotEllipses(kfmeans, kfcovars, pstart,  pmeans, pcovars, pskip)
