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

# This is to add the noise to the measurements! cant add a float to a [float]...
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
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

init_state = [0.5; 450] # initial state
h = 0.01 # time discretisation
tend = 30.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
ys = zeros(N) # only one measurement

f(x, u, w) = Reactor_functions.run_reactor(x, u, h, cstr_model) + w
g(x) = [0.0 1.0]*x# state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 50 #number of particles.
init_state_mean = init_state # initial state mean
init_state_covar = 4.0
init_dist = MvNormal(init_state_mean, init_state_covar) # prior distribution
particles = PF.init_PF(init_dist, nP, 2) # initialise the particles
state_covar = eye(2) # state covariance
state_covar[1] = 1e-4
state_covar[2] = 4.
state_dist = MvNormal(state_covar) # state distribution
meas_covar = eye(1)*4. # measurement covariance
meas_dist = MvNormal(meas_covar) # measurement distribution

fmeans = zeros(2, N)
fcovars = zeros(2,2, N)
# Time step 1
xs[:,1] = init_state
ys[1] = [0.0 1.0]*xs[:, 1] + rand(meas_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys[1], meas_dist, cstr_pf)
fmeans[:,1], fcovars[:,:,1] = PF.getStats(particles)
# Loop through the rest of time
for t=2:N
  xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], 0.0, h, cstr_model) # actual plant
  ys[t] = [0.0 1.0]*xs[:, t] + rand(meas_dist) # measured from actual plant
  PF.filter!(particles, 0.0, ys[t], state_dist, meas_dist, cstr_pf)
  fmeans[:,t], fcovars[:,:,t] = PF.getStats(particles)
end

skip = 20
figure(1) # Kalman Filter Demonstration
x1, = plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
f1, = plot(fmeans[1, 1:skip:end][:], fmeans[2, 1:skip:end][:], "rx", markersize=5, markeredgewidth = 2)
b1 = 0.0
for k=1:skip:N
  p1, p2 = Confidence.plot95(fmeans[:,k], fcovars[:,:, k])
  b1, = plot(p1, p2, "b")
end
ylabel("Temperature [K]")
xlabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1,f1, b1],["Nonlinear Model","Particle Filter Mean", L"Particle Filter $1\sigma$-Ellipse"], loc="best")

figure(2) # Plot filtered results
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
k1, = plot(ts, fmeans[1,:]', "r--", linewidth=3)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1, k1],["Nonlinear Model","Filtered Mean"], loc="best")
xlim([0, tend])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
y2, = plot(ts[1:10:end], ys[1:10:end], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts, fmeans[2,:]', "r--", linewidth=3)
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y2],["Nonlinear Model Measured"], loc="best")
xlim([0, tend])
rc("font",size=22)
