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

# Add a definition for convert to make our lives easier!
# But be careful now!
function Base.convert(::Type{Float64}, x::Array{Float64, 1})
  return x[1]
end

# Specify the nonlinear model
cstr1 = begin
  V = 0.1 #m3
  R = 8.314 #kJ/kmol.K
  CA0 = 1.0 #kmol/m3
  TA0 = 310.0 #K
  dH = -4.78e4 #kJ/kmol
  k0 = 72.0e9 #1/min
  E = 8.314e4 #kJ/kmol
  Cp = 0.239 #kJ/kgK
  rho = 1000.0 #kg/m3
  F = 100e-3 #m3/min
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

cstr2 = begin # slower reaction rate
  V = 0.1 #m3
  R = 8.314 #kJ/kmol.K
  CA0 = 1.0 #kmol/m3
  TA0 = 310.0 #K
  dH = -4.78e4 #kJ/kmol
  k0 = 72.0e9*0.7 #1/min # changed here!
  E = 8.314e4 #kJ/kmol
  Cp = 0.239*1.2 #kJ/kgK # changed here!
  rho = 1000.0 #kg/m3
  F = 100e-3 #m3/min
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

initial_states = [0.57; 410] # initial state
h = 0.001 # time discretisation
tend = 3.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
ys = zeros(2, N)

newC = eye(2)
R = eye(2)
R[1] = 1e-4
R[4] = 2.0
Q = eye(2)
Q[1] = 1e-6
Q[4] = 2.0

A = [0.9 0.1;0.1 0.9]
fun1(x,u,w) = Reactor_functions.run_reactor(x, u, h, cstr1)
fun2(x,u,w) = Reactor_functions.run_reactor(x, u, h, cstr2)
gs(x) = newC*x
F = [fun1, fun2]
G = [gs, gs]

ydists = [MvNormal(R);MvNormal(R)]
xdists = [MvNormal(Q); MvNormal(Q)]
cstr_filter = SPF.Model(F, G, A, xdists, ydists)

nP = 1000
initial_covar = eye(2)
initial_covar[1] = 1e-6
initial_covar[4] = 1.0
xdist = MvNormal(initial_states, initial_covar)
sdist = Categorical([0.5, 0.5])
particles = SPF.init_SPF(xdist, sdist, nP, 2)
fmeans = zeros(2, N)
fcovars = zeros(2,2, N)

us = zeros(N)
measurements = MvNormal(R)
xs[:,1] = initial_states
ys[:, 1] = newC*xs[:, 1] + rand(measurements) # measured from actual plant
SPF.init_filter!(particles, 0.0, ys[:, 1], cstr_filter)
fmeans[:,1], fcovars[:,:,1] = SPF.getStats(particles)
# Loop through the rest of time
for t=2:N
  if ts[t] < 0.1
    xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], us[t-1], h, cstr1) # actual plant
  else
    xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], us[t-1], h, cstr2)
  end
  ys[:, t] = newC*xs[:, t] + rand(measurements) # measured from actual plant
  SPF.filter!(particles, us[t-1], ys[:, t], cstr_filter)
  fmeans[:,t], fcovars[:,:,t] = SPF.getStats(particles)
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
