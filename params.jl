# All the common model parameters go here

# All the required modules
using LLDS
using Reactor
using Ellipse
using LQR
using PF
using PSO
using RBPF
using SPF
using Results
using NLsolve

# Extend the Base Library
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

# Specify the nonlinear model of the broken plant
cstr_model_broken = begin
  V = 5.0 #m3
  R = 8.314 #kJ/kmol.K
  CA0 = 1.0 #kmol/m3
  TA0 = 310.0 #K
  dH = -4.78e4 #kJ/kmol
  k0 = 0.5*72.0e7 #1/min <= the catalyst has denatured
  E = 8.314e4 #kJ/kmol
  Cp = 0.239 #kJ/kgK
  rho = 1000.0 #kg/m3
  F = 100e-3 #m3/min
  Reactor.reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

# Discretise the system
h = 0.1 # time discretisation
tend = 50.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N) # nonlinear plant
linxs = zeros(2, N) # linear plant
xsnofix = zeros(2, N) # broken plant
ys1 = zeros(N) # only measure temperature
ys2 = zeros(2, N) # measure both concentration and temperature
us = zeros(N) # controller input

init_state_covar = eye(2) # prior covariance
init_state_covar[1] = 1e-3
init_state_covar[4] = 4.

pfmeans = zeros(2, N) # Particle Filter means
pfcovars = zeros(2,2, N) # Particle Filter covariances (assumed Gaussian)
rbpfmeans = zeros(2, N) # RBPF means
rbpfcovars = zeros(2,2, N) # RBPF covariances
kfmeans = zeros(2, N) # Kalman Filter means
kfcovars = zeros(2,2, N) # Kalman Filter covariances


# Noise settings
Q = eye(2) # plant noise
Q[1] = 1e-06
Q[4] = 0.1

R1 = eye(1)*10.0 # measurement noise (only temperature)
R2 = eye(2) # measurement noise (both concentration and temperature)
R2[1] = 1e-3
R2[4] = 10.0

# Measurement settings
C2 = eye(2) # we measure both concentration and temperature
C1 = [0.0 1.0] # we measure only temperature

# Controller settings (using quadratic cost function)
QQ = zeros(2, 2)
QQ[1] = 10000.0
RR = 0.00001
