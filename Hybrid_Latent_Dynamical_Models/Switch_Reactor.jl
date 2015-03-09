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

# Specify the first linear model
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
  C = zeros(1,2)
  C[2] = 1.0 #measure temperature
  Q = eye(2) # plant mismatch/noise
  Q[1] = 1e-6
  Q[4] = 4.
  R = 2.0 # measurement noise
  LLDS_functions.LLDS{Float64}(A, B, b, C, Q, R)
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
  C = zeros(1,2)
  C[2] = 1.0 # measure temperature
  Q = eye(2) # plant mismatch/noise
  Q[1] = 1e-6
  Q[4] = 4.
  R = 2.0 # measurement noise
  LLDS_functions.LLDS{Float64}(A, B, b, C, Q, R)
end
fun3(x, u, w) = lin3.A*x + lin3.B*u + lin3.b + w
gun3(x) = lin3.C*x

A = [0.8 0.3 0.01;0.19 0.4 0.19;0.01 0.3 0.8]
switch = SPF.Switch(A)
F = [fun1, fun2, fun3]
G = [gun1, gun2, gun3]
ydists = [MvNormal(eye(1)*lin1.R);MvNormal(eye(1)*lin2.R);MvNormal(eye(1)*lin3.R)]
cstr = SPF.Model(F,G, ydists)

nP = 50
initial_states = [0.57, 395]
initial_covar = eye(2)
initial_covar[1] = 1e-6
initial_covar[4] = 1.0
xdist = MvNormal(initial_states, initial_covar)
sdist = Categorical([0.1, 0.8, 0.1])
particles = SPF.init_SPF(xdist, sdist, nP, 2)
