# Test the Particle Filter.
# We conduct the tests by comparing the posterior
# filtered densities to the analytic Kalman Filter
# solution as calculated by the functions in the
# Linear_Latent_Dynamical_Models folder.

using PyPlot
using Distributions
import PF
cd("..\\CSTR_Model")
using Reactor_functions
cd("..\\Linear_Latent_Dynamical_Models")
using Confidence
cd("..\\Nonlinear_Latent_Dynamical_Models")

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
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

init_state = [0.5; 400] # initial state
h = 0.1 # time discretisation
tend = 150.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
ys = zeros(2, N) # only one measurement

C = eye(2)
f(x, u, w) = Reactor_functions.run_reactor(x, u, h, cstr_model) + w
g(x) = C*x # state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 100 #number of particles.
init_state_mean = init_state # initial state mean
init_state_covar = eye(2)*1e-3 # initial covariance
init_state_covar[4] = 4.0
init_dist = MvNormal(init_state_mean, init_state_covar) # prior distribution
particles = PF.init_PF(init_dist, nP, 2) # initialise the particles
state_covar = eye(2) # state covariance
state_covar[1] = 1e-5
state_covar[2] = 4.
state_dist = MvNormal(state_covar) # state distribution
meas_covar = eye(2)
meas_covar[1] = 1e-3
meas_covar[4] = 10.
meas_dist = MvNormal(meas_covar) # measurement distribution

fmeans = zeros(2, N)
fcovars = zeros(2,2, N)
# Time step 1
xs[:,1] = init_state
ys[:, 1] = C*xs[:, 1] + rand(meas_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys[:, 1], meas_dist, cstr_pf)
fmeans[:,1], fcovars[:,:,1] = PF.getStats(particles)
# Loop through the rest of time
for t=2:N
  xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], 0.0, h, cstr_model) + rand(state_dist) # actual plant
  ys[:, t] = C*xs[:, t] + rand(meas_dist) # measured from actual plant
  PF.filter!(particles, 0.0, ys[:, t], state_dist, meas_dist, cstr_pf)
  fmeans[:,t], fcovars[:,:,t] = PF.getStats(particles)
end
rc("font", family="serif", size=24)
skip = 150
figure(1)
x1, = plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
f1, = plot(fmeans[1, 1:skip:end][:], fmeans[2, 1:skip:end][:], "rx", markersize=5, markeredgewidth = 2)
b1 = 0.0
for k=1:skip:N
  p1, p2 = Confidence.plot95(fmeans[:,k], fcovars[:,:, k])
  b1, = plot(p1, p2, "b")
end
plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
ylabel("Temperature [K]")
xlabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1,f1, b1],["Nonlinear Model","Particle Filter Mean", L"Particle Filter $1\sigma$-Ellipse"], loc="best")

skipm = 20
figure(2) # Plot filtered results
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
y1, = plot(ts[1:skipm:end], ys[1, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k1, = plot(ts, fmeans[1,:]', "r--", linewidth=3)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1, k1],["Nonlinear Model","Filtered Mean"], loc="best")
xlim([0, tend])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
y2, = plot(ts[1:skipm:end], ys[2, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts, fmeans[2,:]', "r--", linewidth=3)
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y2],["Nonlinear Model Measured"], loc="best")
xlim([0, tend])
