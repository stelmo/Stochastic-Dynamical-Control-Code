# Controller using the linear reactor model measuring both concentration and temperature.

tend = 40
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
H = [1.0 0.0] # only attempt to control the concentration
x_off, usp = LQR.offset(A,B,C2,H, ysp) # control offset

f(x, u, w) = A*x + B*u + w
g(x) = C2*x # state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 5000 # number of particles.
prior_dist = MvNormal(init_state-b, init_state_covar) # prior distribution
particles = PF.init_PF(prior_dist, nP, 2) # initialise the particles
state_noise_dist = MvNormal(Q) # state distribution
meas_noise_dist = MvNormal(R2) # measurement distribution

# First time step of the simulation
xs[:,1] = init_state - b # set simulation starting point to the random initial state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measure from actual plant
PF.init_filter!(particles, 0.0, ys2[:, 1], meas_noise_dist, cstr_pf)
pfmeans[:,1], pfcovars[:,:,1] = PF.getStats(particles)

# Setup MPC
horizon = 150
# add state constraints
aline = 10. # slope of constraint line ax + by + c = 0
cline = -411.0 # negative of the y axis intercept
bline = 1.0

us[1] = MPC.mpc_var(pfmeans[:, 1], pfcovars[:,:, 1], horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 10000.0, 1000.0, false, 1.0, Q, 13.8155, true)# get the controller input
temp_states = zeros(2, nP)
tic()
for t=2:N
  xs[:, t] = A*xs[:, t-1] + B*us[t-1] + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measure from actual plant
  PF.filter!(particles, us[t-1], ys2[:, t], state_noise_dist, meas_noise_dist, cstr_pf)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles)

  if t%10 == 0
    us[t] = MPC.mpc_var(pfmeans[:, t], pfcovars[:, :, t], horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 10000.0, 1000.0, false, 1.0, Q, 13.8155, true)
  else
    us[t] = us[t-1]
  end
  if ts[t] in [0.0:2.0:tend]
    Auxiliary.showEstimatedDensity(particles.x .+b, particles.w, temp_states)
  end

end
toc()
pfmeans = pfmeans .+ b
xs = xs .+ b
ys2 = ys2 .+ b


# # Plot the results
Results.plotTracking(ts, xs, ys2, pfmeans, us, 2, ysp+b[1])

Results.plotEllipses(ts, xs, pfmeans, pfcovars, "MPC", [aline, cline], linsystems[2].op, true, 4.6052, 1, "upper left")
