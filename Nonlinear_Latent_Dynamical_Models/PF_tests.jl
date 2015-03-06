# Test the Particle Filter.
# We conduct the tests by comparing the posterior
# filtered densities to the analytic Kalman Filter
# solution as calculated by the functions in the
# Linear_Latent_Dynamical_Models folder.

using PyPlot
import PF
reload("PF.jl")
cd("..\\CSTR_Model")
using Reactor_functions
cd("..\\Nonlinear_Latent_Dynamical_Models")

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
tend = 1.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
xs[:,1] = init_state
ys = zeros(N) # only one measurement

f(x, u, w) = Reactor_functions.run_reactor(x, u, h, cstr_model) + w
g(x) = [0.0 1.0]*x# state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 5 #number of particles.
init_state_mean = init_state # initial state mean
init_state_covar = eye(2)*1e-8 # initial covariance
init_dist = MvNormal(init_state_mean, init_state_covar) # prior distribution
particles = PF.init_PF(init_dist, nP, 1) # initialise the particles

# PF Algorithm
state_covar = eye(2)*1e-6 # initial covariance
state_dist = MvNormal(state_covar) # state distribution
meas_covar = eye(2)*1e-6
meas_dist = MvNormal(meas_covar) # measurement distribution
