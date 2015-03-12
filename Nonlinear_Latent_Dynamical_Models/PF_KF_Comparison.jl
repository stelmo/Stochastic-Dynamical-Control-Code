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
using LLDS_functions
cd("..\\Nonlinear_Latent_Dynamical_Models")

# Specify the nonlinear model
cstr = begin
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

init_state = [0.50; 400]
h = 0.001 # time discretisation
tend = 20.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
xs[:,1] = init_state
ys = zeros(2, N) # only one measurement

xspace = [0.0, 1.0]
yspace = [250, 550]

# Specify the linear model
npoints = 1
linsystems = Reactor_functions.getLinearSystems_randomly(npoints, xspace, yspace, h, cstr) # doesnt work weirdly...
A = linsystems[3].A
B = linsystems[3].B
b = linsystems[3].b
C = eye(2)
Q = eye(2) # plant mismatch/noise
Q[1] = 1e-6
Q[4] = 4.
R = eye(2)
R[1] = 1e-4
R[4] = 10.0 # measurement noise
lin_cstr = LLDS_functions.LLDS(A, B, b, C, Q, R)

f(x, u, w) = A*x + B*u + b + w
g(x) = C*x# state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 500 #number of particles.
init_state_mean = init_state # initial state mean
init_state_covar = eye(2)*1e-6 # initial covariance
init_state_covar[4] = 2.0
init_dist = MvNormal(init_state_mean, init_state_covar) # prior distribution
particles = PF.init_PF(init_dist, nP, 2) # initialise the particles
state_covar = eye(2) # state covariance
state_covar[1] = 1e-4
state_covar[2] = 4.
state_dist = MvNormal(state_covar) # state distribution
meas_covar = eye(2)
meas_covar[1] = 1e-4
meas_covar[4] = 10.
meas_dist = MvNormal(meas_covar) # measurement distribution

fmeans = zeros(2, N)
fcovars = zeros(2,2, N)
filtermeans = zeros(2, N)
filtercovars = zeros(2, 2, N)
# Time step 1
xs[:,1] = init_state
ys[:, 1] = C*xs[:, 1] + rand(meas_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys[:, 1], meas_dist, cstr_pf)
fmeans[:,1], fcovars[:,:,1] = PF.getStats(particles)
filtermeans[:, 1], filtercovars[:, :, 1] = LLDS_functions.init_filter(init_state_mean, init_state_covar, ys[:, 1], lin_cstr)
# Loop through the rest of time
for t=2:N
  xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], 0.0, h, cstr_model) # actual plant
  ys[:, t] = C*xs[:, t] + rand(meas_dist) # measured from actual plant
  PF.filter!(particles, 0.0, ys[:, t], state_dist, meas_dist, cstr_pf)
  fmeans[:,t], fcovars[:,:,t] = PF.getStats(particles)
  filtermeans[:, t], filtercovars[:,:, t] = LLDS_functions.step_filter(filtermeans[:, t-1], filtercovars[:,:, t-1], 0.0, ys[:, t], lin_cstr)
end

skip = 50
figure(1) #  Filter Demonstration
x1, = plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
f1, = plot(fmeans[1, 1:skip:end][:], fmeans[2, 1:skip:end][:], "rx", markersize=5, markeredgewidth = 2)
f2, = plot(filtermeans[1, 1:skip:end][:], filtermeans[2, 1:skip:end][:], "mx", markersize=5, markeredgewidth = 2)
b1 = 0.0
b2 = 0.0
for k=1:skip:N
  p1, p2 = Confidence.plot95(fmeans[:,k], fcovars[:,:, k])
  b1, = plot(p1, p2, "b")

  p3, p4 = Confidence.plot95(filtermeans[:,k], filtercovars[:,:, k])
  b2, = plot(p3, p4, "g")
end
ylabel("Temperature [K]")
xlabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1,f1,f2, b1, b2],["Nonlinear Model","Particle Filter Mean","Kalman Filter Mean", L"Particle Filter $1\sigma$-Ellipse", L"Kalman Filter $1\sigma$-Ellipse"], loc="best")

skipm = 10
skip = 50
figure(2) # Plot filtered results
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
k1, = plot(ts[1:skip:end], fmeans[1,1:skip:end]', "r--", linewidth=3)
y2, = plot(ts[1:skipm:end], ys[1, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k12, = plot(ts[1:skip:end], filtermeans[1, 1:skip:end]', "mo")
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1, k1],["Nonlinear Model","Particle Filter"], loc="best")
xlim([0, tend])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
y2, = plot(ts[1:skipm:end], ys[2, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts[1:skip:end], fmeans[2,1:skip:end]', "r--", linewidth=3)
k22, = plot(ts[1:skip:end], filtermeans[2, 1:skip:end]', "mo")
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y2, k22],["Nonlinear Model Measured", "Kalman Filter"], loc="best")
xlim([0, tend])
rc("font",size=22)
