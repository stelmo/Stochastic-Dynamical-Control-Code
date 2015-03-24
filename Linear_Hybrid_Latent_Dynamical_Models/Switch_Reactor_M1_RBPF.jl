# Implement the augmented switching dynamical system
using PyPlot
using Distributions

import RBPF
reload("RBPF.jl")
import SPF
reload("SPF.jl")
cd("..\\CSTR_Model")
using Reactor_functions
cd("..\\Linear_Latent_Dynamical_Models")
using Confidence
cd("..\\Linear_Hybrid_Latent_Dynamical_Models")

# Add a definition for convert to make our lives easier!
# But be careful now!
function Base.convert(::Type{Float64}, x::Array{Float64, 1})
  return x[1]
end

# Specify the nonlinear model
cstr_model = begin
  V = 5.0; #m3
  R = 8.314; #kJ/kmol.K
  CA0 = 1.0; #kmol/m3
  TA0 = 310.0; #K
  dH = -4.78e4; #kJ/kmol
  k0 = 72.0e7; #1/min
  E = 8.314e4; #kJ/kmol
  Cp = 0.239; #kJ/kgK
  rho = 1000.0; #kg/m3
  F = 100e-3; #m3/min
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

h = 0.001 # time discretisation
tend = 3. # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
ys = zeros(N) # only one measurement

init_state = [0.6; 450] # initial state
C =  [0.0 1.0]
R = eye(1)
R[1] = 4.0
Q = eye(2)
Q[1] = 1e-5
Q[4] = 4.0

# Divide state space into sectors: n by m
nX = 5 # rows
nY = 5 # cols
xspace = [0.0, 1.0]
yspace = [250, 650]

linsystems = Reactor_functions.getLinearSystems(nX, nY, xspace, yspace, h, cstr_model)
# linsystems = Reactor_functions.getLinearSystems_randomly(npoints, xspace, yspace, h, cstr_model)
models, A = RBPF.setup_RBPF(linsystems, C, Q, R)
nP = 500
initial_states = init_state
initial_covar = eye(2)
initial_covar[1] = 1e-6
initial_covar[4] = 1.0
sguess =  SPF.getInitialSwitches(initial_states, linsystems)
particles = RBPF.init_RBPF(Categorical(sguess), initial_states, initial_covar, 2, nP)

fmeans = zeros(2, N)
#
us = zeros(N)
measurements = MvNormal(R)
xs[:,1] = initial_states
ys[1] = C*xs[:, 1] + rand(measurements) # measured from actual plant
RBPF.init_filter!(particles, 0.0, ys[1], models)
fmeans[:,1] = RBPF.getStats(particles)
# Loop through the rest of time
for t=2:N
  xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) # actual plant
  ys[t] = C*xs[:, t] + rand(measurements) # measured from actual plant
  RBPF.filter!(particles, us[t-1], ys[t], models, A)
  fmeans[:, t] = RBPF.getStats(particles)

end


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
