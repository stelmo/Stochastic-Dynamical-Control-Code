# Test the Particle Filter.

using Distributions
using Base.Test
import PF

dircontent = readdir()
if "A.csv" in dircontent
  cd("..\\CSTR_Model")
  import Reactor_functions
  cd("..\\Nonlinear_Latent_Dynamical_Models")
  A = readcsv("A.csv")
  B = readcsv("B.csv")
  b = readcsv("b.csv")
  kfmeans = readcsv("KFmeans.csv")
else
  import Reactor_functions
  A = readcsv(string(pwd(),"/Nonlinear_Latent_Dynamical_Models/A.csv"))
  B = readcsv(string(pwd(),"/Nonlinear_Latent_Dynamical_Models/B.csv"))
  b = readcsv(string(pwd(),"/Nonlinear_Latent_Dynamical_Models/b.csv"))
  kfmeans = readcsv(string(pwd(),"/Nonlinear_Latent_Dynamical_Models/KFmeans.csv"))
end

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

init_state = [0.50; 400]
h = 0.01 # time discretisation
tend = 20.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
xs[:,1] = init_state
ys = zeros(2, N) # only one measurement


C = eye(2)
Q = eye(2) # plant mismatch/noise
Q[1] = 1e-6
Q[4] = 4.
R = eye(2)
R[1] = 1e-4
R[4] = 10.0 # measurement noise

f(x, u, w) = A*x + B*u + b + w
g(x) = C*x# state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 50 #number of particles.
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
# Time step 1
xs[:,1] = init_state
ys[:, 1] = C*xs[:, 1] + rand(meas_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys[:, 1], meas_dist, cstr_pf)
fmeans[:,1], fcovars[:,:,1] = PF.getStats(particles)
# Loop through the rest of time
for t=2:N
  xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], 0.0, h, cstr_model) # actual plant
  ys[:, t] = C*xs[:, t] + rand(meas_dist) # measured from actual plant
  PF.filter!(particles, 0.0, ys[:, t], state_dist, meas_dist, cstr_pf)
  fmeans[:,t], fcovars[:,:,t] = PF.getStats(particles)
end

# Run the tests
tol = 10.0
pf_handler(r::Test.Success) = println("Successful Particle Filter test!")
pf_handler(r::Test.Failure) = error("Failure with the Particle Filter test: $(r.expr)")
pf_handler(r::Test.Error) = rethrow(r)
Test.with_handler(pf_handler) do
  @test maximum(abs(fmeans-kfmeans)) < tol
end
