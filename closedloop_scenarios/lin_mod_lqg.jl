# Controller using the linear reactor model measuring both concentration and temperature.

tend = 100
include("closedloop_params.jl") # load all the parameters and modules

# Get the linear model
linsystems = Reactor.getNominalLinearSystems(h, cstr_model) # cstr_model comes from params.jl
opoint = 2 # the specific operating point we are going to use for control

init_state = [0.55, 450] # random initial point near operating point

# Set the state space model
A = linsystems[opoint].A
B = linsystems[opoint].B
b = linsystems[opoint].b # offset from the origin

# Set point
ysp = linsystems[2].op[1] - b[1] # Medium concentration

# Create the controller
H = [1.0 0.0] # only attempt to control the concentration
x_off, u_off = LQR.offset(A,B,C2,H, ysp) # control offset
K = LQR.lqr(A, B, QQ, RR) # controller

# Set up the KF
kf_cstr = LLDS.llds(A, B, C2, Q, R2) # set up the KF object (measuring both states)
state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(R2)

# First time step of the simulation
xs[:,1] = init_state # set simulation starting point to the random initial state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measure from actual plant
kfmeans[:, 1], kfcovars[:,:, 1] = LLDS.init_filter(init_state-b, init_state_covar, ys2[:, 1]-b, kf_cstr) # filter

us[1] = -K*(kfmeans[:, 1] - x_off) + u_off # controller action

for t=2:N
  xs[:, t] = A*(xs[:, t-1]-b) + B*us[t-1] + b + rand(state_noise_dist) # actual plant

  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measure from actual plant
  kfmeans[:, t], kfcovars[:,:, t] = LLDS.step_filter(kfmeans[:, t-1], kfcovars[:,:, t-1], us[t-1], ys2[:, t]-b, kf_cstr)

  # Compute controller action
  if t%10 == 0
    us[t] = -K*(kfmeans[:, t] - x_off) + u_off # controller action
  else
    us[t] = us[t-1]
  end
end
kfmeans = kfmeans .+ b

# Plot the results
Results.plotTracking(ts, xs, ys2, kfmeans, us, 2, ysp+b[1])
