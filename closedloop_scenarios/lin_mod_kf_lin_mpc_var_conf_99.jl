# Linear Plant controlled with a linear MPC using a KF to estimate the state.
# Stochastic constraints.

tend = 20 # end time of simulation
include("closedloop_params.jl") # load all the parameters and modules

#Get the linear model
linsystems = Reactor.getNominalLinearSystems(h, cstr_model) # cstr_model comes from params.jl
opoint = 2 # the specific operating point we are going to use for control

init_state = [0.5, 450] # random initial point near operating point

# Set the state space model
A = linsystems[opoint].A
B = linsystems[opoint].B
b = linsystems[opoint].b # offset from the origin

# Set point
ysp = linsystems[2].op[1] - b[1] # Medium concentration

# Set up the KF
kf_cstr = LLDS.llds(A, B, C2, Q, R2) # set up the KF object (measuring both states)
state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(R2)

# First time step of the simulation
xs[:,1] = init_state - b # set simulation starting point to the random initial state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measure from actual plant
kfmeans[:, 1], kfcovars[:,:, 1] = LLDS.init_filter(init_state-b, init_state_covar, ys2[:, 1], kf_cstr) # filter

# Setup MPC
horizon = 150
# add state constraints
aline = 10. # slope of constraint line ax + by + c = 0
cline = -408.0 # negative of the y axis intercept
bline = 1.0

us[1] = MPC.mpc_var(kfmeans[:, 1], kfcovars[:,:, 1], horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, 15000.0, 3000.0, false, 1.0, Q, 9.21) # get the controller input
tic()
for t=2:N
  xs[:, t] = A*xs[:, t-1] +  B*us[t-1] + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measure from actual plant
  kfmeans[:, t], kfcovars[:,:, t] = LLDS.step_filter(kfmeans[:, t-1], kfcovars[:,:, t-1], us[t-1], ys2[:, t], kf_cstr)


  us[t] = MPC.mpc_var(kfmeans[:, t], kfcovars[:,:, t], horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, 15000.0, 3000.0, false, 1.0, Q, 9.21) # get the controller input
end
toc()
kfmeans = kfmeans .+ b
xs = xs .+ b
ys2 = ys2 .+ b

# # Plot the results
Results.plotTracking(ts, xs, ys2, kfmeans, us, 2, ysp+b[1])
Results.plotEllipses(ts, xs, kfmeans, kfcovars, "MPC", [aline, cline], linsystems[2].op, true, 9.21)
