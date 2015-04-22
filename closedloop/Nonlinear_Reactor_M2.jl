# Control using the nonlinear reactor model measuring both concentration and temperature

include("../params.jl") # load all the parameters and modules

# Set up the PF
f(x, u, w) = Reactor.run_reactor(x, u, h, cstr_model) + w # transistion function
g(x) = C2*x # measurement function
pf_cstr = PF.Model(f,g) # PF object

# Simulation setup
init_state = [0.5; 400] # initial state

# Initialise the PF
pfmeans = zeros(2, N)
pfcovars = zeros(2,2, N)
nP = 20 # number of particles.
init_state_mean = init_state # initial state mean
init_state_covar = eye(2)*1e-3 # initial covariance
init_state_covar[4] = 4.0
init_dist = MvNormal(init_state_mean, init_state_covar) # prior distribution
particles = PF.init_PF(init_dist, nP, 2) # initialise the particles
state_covar = eye(2) # state covariance
state_covar[1] = 1e-5
state_covar[2] = 0.1
state_noise_dist = MvNormal(state_covar) # state distribution
meas_covar = eye(2)
meas_covar[1] = 1e-3
meas_covar[4] = 10.
meas_noise_dist = MvNormal(meas_covar) # measurement distribution

# Controller setup
y_ca = 0.01 # concentration set point (0.48934869384882873)
offres = PSO.offset(y_ca, cstr_model)
ys2p = [y_ca, offres.zero[1]]
usp = offres.zero[2]

# PSO setup
skip = 7 # hold the controller output for this long
predictionHorizon = 10 # number of controller moves predicted/optimised (control horizon)
swarmSize = 50 # number of particles in the swarm

# Time step 1
xs[:,1] = init_state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys2[:, 1], meas_noise_dist, pf_cstr)
pfmeans[:,1], pfcovars[:,:,1] = PF.getStats(particles)
# Controlle action
swarm, sol = PSO.initswarm(swarmSize, predictionHorizon, -1000.0, 1000.0, particles, ys2p, usp, QQ, RR, state_noise_dist, cstr_model, skip, h)
us[1] = PSO.optimise!(swarm, sol, particles, ys2p, usp, QQ, RR, state_noise_dist, cstr_model, skip, h)

for t=2:50
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  PF.filter!(particles, us[t-1], ys2[:, t], state_noise_dist, meas_noise_dist, pf_cstr)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles)

  # Controller action
  if t%skip == 0
    swarm, sol = PSO.initswarm(swarmSize, predictionHorizon, -1000.0, 1000.0, particles, ys2p, usp, QQ, RR, state_noise_dist, cstr_model, skip, h)
    us[t] = PSO.optimise!(swarm, sol, particles, ys2p, usp, QQ, RR, state_noise_dist, cstr_model, skip, h)
    # println("Time: ", ts[t-1])
  end
  if t%skip != 0
    us[t] = us[t-1]
  end
  # println("Controller Input: ", us[t-1])

end

#  Plot Results
rc("font", family="serif", size=24)
skipmeas = int(length(ts)/20)

figure(1)
subplot(3,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
y1, = plot(ts[1:skipmeas:end], ys2[1, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
k1, = plot(ts, pfmeans[1,:]', "rx", markersize=5, markeredgewidth=2)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1, k1],["Nonlinear Model","Filtered Mean"], loc="best")
xlim([0, tend])
subplot(3,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
y2, = plot(ts[1:skipmeas:end], ys2[2, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts, pfmeans[2,:]', "rx", markersize=5, markeredgewidth=2)
ylabel("Temperature [K]")
legend([y2],["Nonlinear Model Measured"], loc="best")
xlim([0, tend])
subplot(3,1,3)
plot(ts, us)
xlabel("Time [min]")
ylabel("Controller Input")
