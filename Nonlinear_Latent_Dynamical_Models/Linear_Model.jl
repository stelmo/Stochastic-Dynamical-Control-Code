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

# Add a definition for convert to make our lives easier!
# But be careful now!
function Base.convert(::Type{Float64}, x::Array{Float64, 1})
  return x[1]
end

# Specify the nonlinear model
cstr_model = begin
  V = 0.1; #m3
  R = 8.314; #kJ/kmol.K
  CA0 = 1.0; #kmol/m3
  TA0 = 310.0; #K
  dH = -4.78e4; #kJ/kmol
  k0 = 72.0e9; #1/min
  E = 8.314e4; #kJ/kmol
  Cp = 0.239; #kJ/kgK
  rho = 1000.0; #kg/m3
  F = 100e-3; #m3/min
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

init_state = [0.57; 395] # initial state
h = 0.001 # time discretisation
tend = 2.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
ys = zeros(N) # only one measurement

# Specify the linear model
J = readcsv("J_ss.csv")
ss = readcsv("ss.csv")'[:,1]
lin_cstr = begin
  A = eye(2)+h*J
  B = zeros(2,1)
  B[2] = 1./(cstr_model.rho*cstr_model.Cp*cstr_model.V)
  b = -h*J*ss
  C = zeros(1,2)
  C[2] = 1.0 #measure temperature
  Q = eye(2) # plant mismatch/noise
  Q[1] = 1e-6
  Q[4] = 4.
  R = 2.0 # measurement noise
  LLDS_functions.LLDS{Float64}(A, B, b, C, Q, R)
end

f(x, u, w) = A*x + B*u + b + w
g(x) = C*x# state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 100 #number of particles.
init_state_mean = init_state # initial state mean
init_state_covar = eye(2)*1e-6 # initial covariance
init_state_covar[4] = 2.0
init_dist = MvNormal(init_state_mean, init_state_covar) # prior distribution
particles = PF.init_PF(init_dist, nP, 2) # initialise the particles
state_dist = MvNormal(Q) # state distribution
meas_covar = eye(1)*R # measurement covariance
meas_dist = MvNormal(meas_covar) # measurement distribution

fmeans = zeros(2, N)
fcovars = zeros(2,2, N)
# Time step 1
xs[:,1] = init_state
ys[1] = C*xs[:, 1] + rand(meas_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys[1], state_dist, meas_dist, cstr_pf)
fmeans[:,1], fcovars[:,:,1] = PF.getStats(particles)
filtermeans = zeros(2, N)
filtercovars = zeros(2,2, N)
filtermeans[:, 1], filtercovars[:,:, 1] = LLDS_functions.init_filter(init_state_mean, init_state_covar, ys[1], lin_cstr)
# Loop through the rest of time
for t=2:N
  xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], 0.0, h, cstr_model) # actual plant
  ys[t] = C*xs[:, t] + rand(meas_dist) # measured from actual plant
  PF.filter!(particles, 0.0, ys[t], state_dist, meas_dist, cstr_pf)
  fmeans[:,t], fcovars[:,:,t] = PF.getStats(particles)
  filtermeans[:, t], filtercovars[:,:, t] = LLDS_functions.step_filter(filtermeans[:, t-1], filtercovars[:,:, t-1], 0.0, ys[t], lin_cstr)
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

figure(2) # Plot filtered results
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
k1, = plot(ts[1:skip:end], fmeans[1,1:skip:end]', "r--", linewidth=3)
k12, = plot(ts[1:skip:end], filtermeans[1, 1:skip:end]', "mo")
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1, k1],["Nonlinear Model","Particle Filter"], loc="best")
xlim([0, tend])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
y2, = plot(ts[1:10:end], ys[1:10:end], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts[1:skip:end], fmeans[2,1:skip:end]', "r--", linewidth=3)
k22, = plot(ts[1:skip:end], filtermeans[2, 1:skip:end]', "mo")
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y2, k22],["Nonlinear Model Measured", "Kalman Filter"], loc="best")
xlim([0, tend])
rc("font",size=22)
