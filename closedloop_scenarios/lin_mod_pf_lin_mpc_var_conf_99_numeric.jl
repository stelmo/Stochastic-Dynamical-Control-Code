# Controller using the linear reactor model measuring both concentration and temperature.

tend = 20
include("closedloop_params.jl") # load all the parameters and modules

# Get the linear model
linsystems = Reactor.getNominalLinearSystems(h, cstr_model) # cstr_model comes from params.jl
opoint = 2 # the specific operating point we are going to use for control

init_state = [0.5, 450] # random initial point near operating point

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
nP = 20000 # number of particles.
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
cline = -412.0 # negative of the y axis intercept
bline = 1.0

Ndiv = length([0:3.0:tend])
kldiv = zeros(Ndiv) # Kullback-Leibler Divergence as a function of time
klts = zeros(Ndiv)
ndivcounter = 1
temp_states = zeros(2, nP)

us[1] = MPC.mpc_var(pfmeans[:, 1], pfcovars[:,:, 1], horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 15000.0, 1000.0, false, 1.0, Q, 9.21, true)# get the controller input
kldiv[ndivcounter] = Auxiliary.KL(particles.x, particles.w, pfmeans[:, 1], pfcovars[:,:, 1], temp_states)
klts[ndivcounter] = 0.0
ndivcounter += 1

tic()
for t=2:N
  xs[:, t] = A*xs[:, t-1] + B*us[t-1] + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measure from actual plant
  PF.filter!(particles, us[t-1], ys2[:, t], state_noise_dist, meas_noise_dist, cstr_pf)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles)

  us[t] = MPC.mpc_var(pfmeans[:, t], pfcovars[:, :, t], horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 15000.0, 1000.0, false, 1.0, Q, 9.21, true)

  if ts[t] in [0.0:3.0:tend]
    kldiv[ndivcounter] = Auxiliary.KL(particles.x, particles.w, pfmeans[:, t], pfcovars[:,:, t], temp_states)
    klts[ndivcounter] = ts[t]
    ndivcounter += 1
  end
end
toc()
pfmeans =pfmeans .+ b
xs = xs .+ b
ys2 = ys2 .+ b


# # Plot the results
Results.plotTracking(ts, xs, ys2, pfmeans, us, 2, ysp+b[1])

Results.plotEllipses(ts, xs, pfmeans, pfcovars, "MPC", [aline, cline], linsystems[2].op, true)

Results.plotKLdiv(klts, kldiv)
