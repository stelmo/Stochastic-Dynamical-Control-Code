# Test the Particle Filter.
# We conduct the tests by comparing the posterior
# filtered densities to the analytic Kalman Filter
# solution as calculated by the functions in the
# Linear_Latent_Dynamical_Models folder.

using PF
using Reactor
using PSO


function Base.convert(::Type{Float64}, x::Array{Float64, 1})
  return x[1]
end

# Specify the nonlinear model
cstr_model = begin
  V = 5.0 #m3
  R = 8.314 #kJ/kmol.K
  CA0 = 1.0 #kmol/m3
  TA0 = 310.0 #K
  dH = -4.78e4 #kJ/kmol
  k0 = 72.0e7 #1/min
  E = 8.314e4 #kJ/kmol
  Cp = 0.239 #kJ/kgK
  rho = 1000.0 #kg/m3
  F = 100e-3 #m3/min
  Reactor.reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

h = 0.1 # time discretisation
C = eye(2)
f(x, u, w) = Reactor.run_reactor(x, u, h, cstr_model) + w
g(x) = C*x # state observation
cstr_pf = PF.Model(f,g)


init_state = [0.5; 400] # initial state
tend = 20.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
ys = zeros(2, N) # only one measurement

# Initialise the PF
nP = 20 #number of particles.
init_state_mean = init_state # initial state mean
init_state_covar = eye(2)*1e-3 # initial covariance
init_state_covar[4] = 4.0
init_dist = MvNormal(init_state_mean, init_state_covar) # prior distribution
particles = PF.init_PF(init_dist, nP, 2) # initialise the particles
state_covar = eye(2) # state covariance
state_covar[1] = 1e-5
state_covar[2] = 0.1
state_dist = MvNormal(state_covar) # state distribution
meas_covar = eye(2)
meas_covar[1] = 1e-3
meas_covar[4] = 10.
meas_dist = MvNormal(meas_covar) # measurement distribution

# Control
R = 0.0
Q = [1.0 0.0;0.0 0.0]
y_ca = 0.48934869384882873 # concentration set point
offres = PSO.offset(y_ca, cstr_model)
ysp = [y_ca, offres.zero[1]]
usp = offres.zero[2]

fmeans = zeros(2, N)
fcovars = zeros(2,2, N)
us = zeros(N-1)
# Time step 1
xs[:,1] = init_state
ys[:, 1] = C*xs[:, 1] + rand(meas_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys[:, 1], meas_dist, cstr_pf)

fmeans[:,1], fcovars[:,:,1] = PF.getStats(particles)

skip = 10 # hold the controller output for this long
predictionHorizon = 5 # number of controller moves
swarmSize = 50

# swarm, sol = PSO.initswarm(swarmSize, predictionHorizon, -500.0, 500.0, particles, ysp, usp, Q, R, state_dist, cstr_model, skip, h)
# @time PSO.optimise!(swarm, sol, particles, ysp, usp, Q, R, state_dist, cstr_model, skip, h)

##Loop through the rest of time
tic()
swarm, sol = PSO.initswarm(swarmSize, predictionHorizon, -500.0, 500.0, particles, ysp, usp, Q, R, state_dist, cstr_model, skip, h)
us[1] = PSO.optimise!(swarm, sol, particles, ysp, usp, Q, R, state_dist, cstr_model, skip, h)

for t=2:N
  if t%skip == 0
    swarm, sol = PSO.initswarm(swarmSize, predictionHorizon, -1000.0, 1000.0, particles, ysp, usp, Q, R, state_dist, cstr_model, skip, h)
    us[t-1] = PSO.optimise!(swarm, sol, particles, ysp, usp, Q, R, state_dist, cstr_model, skip, h)
    println("Time: ", ts[t-1])
  end
  if t%skip != 0 && t != 2
    us[t-1] = us[t-2]
  end

  println("Controller Input: ", us[t-1])

  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_dist) # actual plant
  ys[:, t] = C*xs[:, t] + rand(meas_dist) # measured from actual plant
  PF.filter!(particles, us[t-1], ys[:, t], state_dist, meas_dist, cstr_pf)
  fmeans[:,t], fcovars[:,:,t] = PF.getStats(particles)
end
toc()
skipm = 20
rc("font", family="serif", size=24)
skip = 150
figure(1)

subplot(3,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
y1, = plot(ts[1:skipm:end], ys[1, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k1, = plot(ts, fmeans[1,:]', "rx", markersize=5, markeredgewidth=2)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1, k1],["Nonlinear Model","Filtered Mean"], loc="best")
xlim([0, tend])
subplot(3,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
y2, = plot(ts[1:skipm:end], ys[2, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts, fmeans[2,:]', "rx", markersize=5, markeredgewidth=2)
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y2],["Nonlinear Model Measured"], loc="best")
xlim([0, tend])
subplot(3,1,3)
plot(ts[1:end-1], us)
