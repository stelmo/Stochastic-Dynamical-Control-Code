# Nonlinear plant model controlled with a linear MPC. The control goal is to steer
# the system to the unstead operating point. Deterministic contraints.

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

f(x, u, w) = Reactor.run_reactor(x, u, h, cstr_model) + w
g(x) = C2*x # state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 200 # number of particles.
prior_dist = MvNormal(init_state, init_state_covar) # prior distribution
particles = PF.init_PF(prior_dist, nP, 2) # initialise the particles
state_noise_dist = MvNormal(Q) # state distribution
meas_noise_dist = MvNormal(R2) # measurement distribution

# Setup MPC
horizon = 150
# add state constraints
aline = 10. # slope of constraint line ax + by + c = 0
cline = -406.0 # negative of the y axis intercept
bline = 1.0


mcdists = zeros(2, mcN)
xconcen = zeros(N, mcN)
mcerrs = zeros(mcN)
tic()
for mciter=1:mcN
  prior_dist = MvNormal(init_state, init_state_covar) # prior distribution
  particles = PF.init_PF(prior_dist, nP, 2) # initialise the particles

  # First time step of the simulation
  xs[:,1] = init_state # set simulation starting point to the random initial state
  ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measure from actual plant
  PF.init_filter!(particles, 0.0, ys2[:, 1], meas_noise_dist, cstr_pf)
  pfmeans[:,1], pfcovars[:,:,1] = PF.getStats(particles)

  us[1] = MPC.mpc_mean(pfmeans[:, 1]-b, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 20000.0, 1000.0, false)# get the controller input

  for t=2:N
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
    ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measure from actual plant
    PF.filter!(particles, us[t-1], ys2[:, t], state_noise_dist, meas_noise_dist, cstr_pf)
    pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles)
    if t%10 == 0
      us[t] = MPC.mpc_mean(pfmeans[:, t]-b, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 20000.0, 1000.0, false)
    else
      us[t] = us[t-1]
    end
  end
  mcerrs[mciter] = Results.calcError2(xs, ysp+b[1])
  xconcen[:, mciter] = xs[1, :]
  Results.getMCRes!(xs, pfcovars, [aline, cline], mcdists, mciter, h)
end
toc()

# nocount = count(x->x==0.0, mcdists[1,:])
# filteredResults = zeros(2, mcN-nocount)
# counter = 1
# for k=1:mcN
#   if mcdists[1, k] != 0.0
#   filteredResults[:, counter] = mcdists[:, k]
#   counter += 1
#   end
# end
#
# writecsv("nonlinmod_pf_mean.csv", filteredResults)
println("The absolute MC average error is: ", sum(abs(mcerrs))/mcN)
writecsv("nonlinmod_pf_mean_mc2.csv", xconcen)
