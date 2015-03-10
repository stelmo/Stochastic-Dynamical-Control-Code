# Implement the augmented switching dynamical system
using PyPlot
using Distributions
import SPF
reload("SPF.jl")
cd("..\\CSTR_Model")
using Reactor_functions
cd("..\\Linear_Latent_Dynamical_Models")
using Confidence
using LLDS_functions
cd("..\\Hybrid_Latent_Dynamical_Models")

warn("Hardcoded for 2 state system!")

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

init_state = [0.5; 410] # initial state
h = 0.01 # time discretisation
tend = 3.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
ys = zeros(2, N) # only one measurement

newC = eye(2)
R = eye(2)
R[1] = 1e-4
R[4] = 2.0

# Specify the first linear model
J1 = readcsv("J_ss1.csv")
ss1 = readcsv("ss1.csv")'[:,1]
lin1 = begin
  A = eye(2)+h*J1
  B = zeros(2,1)
  B[2] = h/(cstr_model.rho*cstr_model.Cp*cstr_model.V)
  b = -h*J1*ss1
  # C = zeros(1,2)
  # C[2] = 1.0 #measure temperature
  C = newC
  Q = eye(2) # plant mismatch/noise
  Q[1] = 1e-6
  Q[4] = 4.
  # R = 2.0 # measurement noise
  LLDS_functions.LLDS{Array{Float64,2}}(A, B, b, C, Q, R)
end
fun1(x, u, w) = lin1.A*x + lin1.B*u + lin1.b + w
gun1(x) = lin1.C*x

# Specify the second linear model
J2 = readcsv("J_ss2.csv")
ss2 = readcsv("ss2.csv")'[:,1]
lin2 = begin
  A = eye(2)+h*J2
  B = zeros(2,1)
  B[2] = h/(cstr_model.rho*cstr_model.Cp*cstr_model.V)
  b = -h*J2*ss2
  # C = zeros(1,2)
  # C[2] = 1.0 #measure temperature
  C = newC
  Q = eye(2) # plant mismatch/noise
  Q[1] = 1e-6
  Q[4] = 4.
  # R = 2.0 # measurement noise
  LLDS_functions.LLDS{Array{Float64,2}}(A, B, b, C, Q, R)
end
fun2(x, u, w) = lin2.A*x + lin2.B*u + lin2.b + w
gun2(x) = lin2.C*x

# Specify the third linear model
J3 = readcsv("J_ss3.csv")
ss3 = readcsv("ss3.csv")'[:,1]
lin3 = begin
  A = eye(2)+h*J3
  B = zeros(2,1)
  B[2] = h/(cstr_model.rho*cstr_model.Cp*cstr_model.V)
  b = -h*J3*ss3
  # C = zeros(1,2)
  # C[2] = 1.0 # measure temperature
  C = newC
  Q = eye(2) # plant mismatch/noise
  Q[1] = 1e-6
  Q[4] = 4.
  # R = 2.0 # measurement noise
  LLDS_functions.LLDS{Array{Float64,2}}(A, B, b, C, Q, R)
end
fun3(x, u, w) = lin3.A*x + lin3.B*u + lin3.b + w
gun3(x) = lin3.C*x

A = [0.5 0.3 0.1;
     0.4 0.4 0.4;
     0.1 0.3 0.5]
F = [fun1, fun2, fun3]
G = [gun1, gun2, gun3]
ydists = [MvNormal(lin1.R);MvNormal(lin2.R);MvNormal(lin3.R)]
xdists = [MvNormal(lin1.Q); MvNormal(lin2.Q); MvNormal(lin3.Q)]
cstr = SPF.Model(F, G, A, xdists, ydists)

nP = 1000
initial_states = init_state
initial_covar = eye(2)
initial_covar[1] = 1e-6
initial_covar[4] = 1.0
xdist = MvNormal(initial_states, initial_covar)
sdist = Categorical([0.1, 0.8, 0.1])
particles = SPF.init_SPF(xdist, sdist, nP, 2)
fmeans = zeros(2, N)
fcovars = zeros(2,2, N)

us = zeros(N)
measurements = MvNormal(R)
xs[:,1] = initial_states
ys[:, 1] = newC*xs[:, 1] + rand(measurements) # measured from actual plant
SPF.init_filter!(particles, 0.0, ys[:, 1], cstr)
fmeans[:,1], fcovars[:,:,1] = SPF.getStats(particles)
# Loop through the rest of time
for t=2:N
  xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) # actual plant
  ys[:, t] = newC*xs[:, t] + rand(measurements) # measured from actual plant
  SPF.filter!(particles, us[t-1], ys[:, t], cstr)
  fmeans[:,t], fcovars[:,:,t] = SPF.getStats(particles)
  # println("*******")
  # println(count((x)->x==2, particles.s))
  # if ts[t] > 0.2
  #   us[t] = 100.
  # end
end

skip = 50
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
y2, = plot(ts[1:10:end], ys[1, 1:10:end][:], "kx", markersize=5, markeredgewidth=1)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1, k1],["Nonlinear Model","Filtered Mean"], loc="best")
xlim([0, tend])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
y2, = plot(ts[1:10:end], ys[2, 1:10:end][:], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts, fmeans[2,:]', "r--", linewidth=3)
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y2],["Nonlinear Model Measured"], loc="best")
xlim([0, tend])
rc("font",size=22)
