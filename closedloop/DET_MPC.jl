# Controller using the linear reactor model measuring both concentration and temperature.

tend = 80
include("params.jl") # load all the parameters and modules

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
# ysp = linsystems[2].op[1] - b[1] # Medium concentration
ysp = 0.65 - b[1]
H = [1.0 0.0] # only attempt to control the concentration
x_off, usp = LQR.offset(A,B,C2,H, ysp) # control offset


# First time step of the simulation
xs[:,1] = init_state - b # set simulation starting point to the random initial state

# Setup MPC
horizon = 150
# add state constraints
aline = 10. # slope of constraint line ax + by + c = 0
cline = -300.0 # negative of the y axis intercept
bline = 1.0

d = x_off - xs[:, 1]
# res = MPC.mpc_targets(ysp, d, A, B)
usp = 0.0
us[1] = MPC.mpc_mean(xs[:, 1], horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 15000.0, 1000.0, false)# get the controller input
tic()
for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1]+b, us[t-1], h, cstr_model) -b
  # xs[:, t] = A*xs[:, t-1] + B*us[t-1]
  d = xs[:, t] - (A*xs[:, t-1] + B*us[t-1])
  # println(d)
  us[t] = MPC.mpc_mean(xs[:, t], horizon, A, B, b, aline, bline, cline, QQ, RR, ysp+0.002, usp[1], 15000.0, 1000.0, false)
end
toc()
xs = xs .+ b

# Plot the results
Results.plotTracking(ts, xs, xs, xs, us, 2, ysp+b[1])
