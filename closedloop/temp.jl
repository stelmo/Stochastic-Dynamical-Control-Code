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
ysp = linsystems[1].op[1] - b[1] # Low concentration
# ysp = linsystems[2].op[1] - b[1] # Medium concentration
# ysp = linsystems[3].op[1] - b[1] # High concentration
# ysp = 0.1 - b[1]

# Create the controller
H = [1.0 0.0] # only attempt to control the concentration
x_off, u_off = LQR.offset(A,B,C2,H, ysp) # control offset
K = LQR.lqr(A, B, QQ, RR) # controller

# First time step of the simulation
xs[:,1] = init_state - b # set simulation starting point to the random initial state

us[1] = -K*(xs[:, 1] - x_off) + u_off # controller action
d = 0.0
for t=2:N
  # xs[:, t] = Reactor.run_reactor(xs[:, t-1]+b, us[t-1], h, cstr_model)-b
  xs[:, t] = A*xs[:, t-1] + B*us[t-1]

  # Compute controller action
  d = xs[:, t] - (A*xs[:, t-1] + B*us[t-1])
  us[t] = -K*(xs[:, t] - x_off) + u_off # controller action
end
xs = xs .+b
# Plot the results
Results.plotTracking(ts, xs, ys2, xs, us, 2, ysp+b[1])
