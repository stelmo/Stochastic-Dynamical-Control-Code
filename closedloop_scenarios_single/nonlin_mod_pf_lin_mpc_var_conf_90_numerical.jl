# Nonlinear plant model controlled with a linear MPC. The control goal is to steer
# the system to the unstead operating point. Stochastic contraints. Numerical evaluation
# of the Gaussian assumption.

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

f(x, u, w) = Reactor.run_reactor(x, u, h, cstr_model) + w
g(x) = C2*x # state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 5000 # number of particles.
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

Ndiv = length([0:3.0:tend])
kldiv = zeros(Ndiv) # Kullback-Leibler Divergence as a function of time
basediv = zeros(Ndiv) # Baseline
unidiv = zeros(Ndiv) # Uniform comparison
klts = zeros(Ndiv)
ndivcounter = 1
temp_states = zeros(2, nP)

us[1] = MPC.mpc_var(pfmeans[:, 1]-b, pfcovars[:,:, 1], horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 20000.0, 1000.0, false, 1.0, Q, 9.21, true)# get the controller input

kldiv[ndivcounter] = Auxiliary.KL(particles.x, particles.w, pfmeans[:, 1], pfcovars[:,:, 1], temp_states)
basediv[ndivcounter] = Auxiliary.KLbase(pfmeans[:, 1], pfcovars[:,:, 1], temp_states, nP)
unidiv[ndivcounter] = Auxiliary.KLuniform(pfmeans[:, 1], pfcovars[:,:, 1], temp_states, nP)
klts[ndivcounter] = 0.0
ndivcounter += 1

tic()
for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measure from actual plant
  PF.filter!(particles, us[t-1], ys2[:, t], state_noise_dist, meas_noise_dist, cstr_pf)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles)

  if t%10 == 0
    us[t] = MPC.mpc_var(pfmeans[:, t]-b, pfcovars[:, :, t], horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp[1], 20000.0, 1000.0, false, 1.0, Q, 9.21, true)
  else
    us[t] = us[t-1]
  end

  if ts[t] in [0.0:3.0:tend]
    kldiv[ndivcounter] = Auxiliary.KL(particles.x, particles.w, pfmeans[:, t], pfcovars[:,:, t], temp_states)
    basediv[ndivcounter] = Auxiliary.KLbase(pfmeans[:, t], pfcovars[:,:, t], temp_states, nP)
    unidiv[ndivcounter] = Auxiliary.KLuniform(pfmeans[:, t]-b, pfcovars[:,:, t], temp_states, nP)
    klts[ndivcounter] = ts[t]
    ndivcounter += 1
  end
end
toc()

# # Plot the results
Results.plotTracking(ts, xs, ys2, pfmeans, us, 2, ysp+b[1])

Results.plotEllipses(ts, xs, pfmeans, pfcovars, "MPC", [aline, cline], linsystems[2].op, true, 4.6052, 1, "upper right")

Results.plotKLdiv(klts, kldiv, basediv, unidiv, false)
Results.plotKLdiv(klts, log(kldiv), log(basediv), log(unidiv), true)
println("The average divergence for the baseline is: ",1.0/length(klts)*sum(basediv))
println("The average divergence for the approximation is: ",1.0/length(klts)*sum(kldiv))
println("The average divergence for the uniform is: ",1.0/length(klts)*sum(unidiv))
