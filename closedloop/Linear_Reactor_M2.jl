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

# Set point
# ysp = linsystems[1].op[1] - b[1] # Low concentration
ysp = linsystems[2].op[1] - b[1] # Medium concentration
# ysp = linsystems[3].op[1] - b[1] # High concentration
# ysp = 0.1 - b[1]

# Create the controller
H = [1.0 0.0] # only attempt to control the concentration
x_off, u_off = LQR.offset(A,B,C2,H, ysp) # control offset
K = LQR.lqr(A, B, QQ, RR) # controller

# Set up the KF
kf_cstr = LLDS.llds(A, B, C2, Q, R2) # set up the KF object (measuring both states)
state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(R2)
kfmeans = zeros(2, N)
kfcovars = zeros(2,2, N)

# Filter initialisation
init_mean = init_state - b # state space offset
init_covar = eye(2)
init_covar[1] = 1e-3
init_covar[4] = 4.

# First time step of the simulation
xs[:,1] = init_state # set simulation starting point to the random initial state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measure from actual plant
kfmeans[:, 1], kfcovars[:,:, 1] = LLDS.init_filter(init_mean, init_covar, ys2[:, 1]-b, kf_cstr) # filter
us[1] = -K*(kfmeans[:, 1] - x_off) + u_off # controller action

for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measure from actual plant
  kfmeans[:, t], kfcovars[:,:, t] = LLDS.step_filter(kfmeans[:, t-1], kfcovars[:,:, t-1], us[t-1], ys2[:, t]-b, kf_cstr)

  # Compute controller action
  us[t] = -K*(kfmeans[:, t] - x_off) + u_off # controller action

end
kfmeans = kfmeans .+ b

# Plot the results
Results.plotTracking(ts, xs, ys2, kfmeans, us, 2)
