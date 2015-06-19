# Controller using the linear reactor model measuring both concentration and temperature.

tend = 80
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

# Set point
H = [1.0 0.0] # only attempt to control the concentration
x_off, usp = LQR.offset(A,B,C2,H, ysp) # control offset

# Set up the KF
kf_cstr = LLDS.llds(A, B, C2, Q, R2) # set up the KF object (measuring both states)
state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(R2)

# First time step of the simulation
xs[:,1] = init_state # set simulation starting point to the random initial state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measure from actual plant
kfmeans[:, 1], kfcovars[:,:, 1] = LLDS.init_filter(init_state-b, init_state_covar, ys2[:, 1]-b, kf_cstr) # filter

# Setup MPC
horizon = 150
# add state constraints
aline = 10. # slope of constraint line ax + by + c = 0
cline = -406.0 # negative of the y axis intercept
bline = 1.0

us[1] = MPC.mpc_var(kfmeans[:, 1], kfcovars[:,:, 1], horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 20000.0, 1000.0, false, 1.0, Q, 9.21, true) # get the controller input
tic()
for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measure from actual plant
  kfmeans[:, t], kfcovars[:,:, t] = LLDS.step_filter(kfmeans[:, t-1], kfcovars[:,:, t-1], us[t-1], ys2[:, t]-b, kf_cstr)

  if t%10 == 0
    us[t] = MPC.mpc_var(kfmeans[:, t], kfcovars[:,:, t], horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 20000.0, 1000.0, false, 1.0, Q, 9.21, true) # get the controller input
  else
    us[t] = us[t-1]
  end
end
toc()
kfmeans = kfmeans .+ b

# # Plot the results
Results.plotTracking(ts, xs, ys2, kfmeans, us, 2, ysp+b[1])
Results.plotEllipses(ts, xs, kfmeans, kfcovars, "MPC", [aline, cline], linsystems[2].op, true, 4.6052, 1, "upper right")
Results.checkConstraint(ts, xs, [aline, cline])
Results.calcError(xs, ysp+b[1])
Results.calcEnergy(us, 0.0, h)
