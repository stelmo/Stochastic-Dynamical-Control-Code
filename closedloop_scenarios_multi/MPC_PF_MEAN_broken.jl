# Nonlinear plant model controlled with a linear MPC. The control goal is to steer
# the system to the unstead operating point. Deterministic contraints.

tend = 300
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

f(x, u, w) = Reactor.run_reactor(x, u, h, cstr_model) + w
g(x) = C2*x # state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 200 # number of particles.
prior_dist = MvNormal(init_state, init_state_covar) # prior distribution
particles = PF.init_PF(prior_dist, nP, 2) # initialise the particles
state_noise_dist = MvNormal(Q) # state distribution
meas_noise_dist = MvNormal(R2) # measurement distribution

# First time step of the simulation
xs[:,1] = init_state # set simulation starting point to the random initial state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measure from actual plant
PF.init_filter!(particles, 0.0, ys2[:, 1], meas_noise_dist, cstr_pf)
pfmeans[:,1], pfcovars[:,:,1] = PF.getStats(particles)

# Setup MPC
horizon = 150
# add state constraints
aline = 10. # slope of constraint line ax + by + c = 0
cline = -400.0 # negative of the y axis intercept
bline = 1.0

us[1] = MPC.mpc_mean(pfmeans[:, 1]-b, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 15000.0, 1000.0, false)# get the controller input
tic()
d = zeros(2)
for t=2:N
  d = zeros(2)
  if ts[t] < 100
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
    xtemp = xs[:, t-1] - b
    d = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) - (A*xtemp + B*us[t-1] + b)
  else
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model_broken) + rand(state_noise_dist) # actual plant
    xtemp = xs[:, t-1] - b
    d = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model_broken) - (A*xtemp + B*us[t-1] + b)
  end

  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measure from actual plant
  PF.filter!(particles, us[t-1], ys2[:, t], state_noise_dist, meas_noise_dist, cstr_pf)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles)
  if t%10 == 0
    us[t] = MPC.mpc_mean(pfmeans[:, t]-b, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 15000.0, 1000.0, false, d)
  else
    us[t] = us[t-1]
  end
end
toc()

# # Plot the results
Results.plotTracking(ts, xs, ys2, pfmeans, us, 2, ysp+b[1])
Results.plotEllipses(ts, xs, pfmeans, pfcovars, "MPC", [aline, cline], [linsystems[2].op[1], 422.6], true, 4.6052, 1, "upper right")
Results.checkConstraint(ts, xs, [aline, cline])
Results.calcError(xs, ysp+b[1])
Results.calcEnergy(us, 0.0, h)
