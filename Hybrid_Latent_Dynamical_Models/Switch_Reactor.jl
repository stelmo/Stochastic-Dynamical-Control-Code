# Implement the augmented switching dynamical system
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
J1 = readcsv("J_ss1.csv")
ss1 = readcsv("ss1.csv")'[:,1]
lin1 = begin
  A = eye(2)+h*J1
  B = zeros(2,1)
  B[2] = h/(cstr_model.rho*cstr_model.Cp*cstr_model.V)
  b = -h*J1*ss1
  C = zeros(1,2)
  C[2] = 1.0 #measure temperature
  Q = eye(2) # plant mismatch/noise
  Q[1] = 1e-6
  Q[4] = 4.
  R = 2.0 # measurement noise
  LLDS_functions.LLDS{Float64}(A, B, b, C, Q, R)
end
