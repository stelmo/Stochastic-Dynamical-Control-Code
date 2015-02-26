# This script will run the linear Du reactor.

using PyPlot
using LLDS_functions
using Distributions
import DuReactor_functions
reload("DuReactor_functions.jl")

# Specify the system parameters
du = begin
  phi = 0.072
  q = 1.0
  beta = 8.0
  delta = 0.3
  lambda = 20.0
  x1f = 1.0
  x2f = 0.0
  DuReactor_functions.DuReactor(phi, q, beta, delta, lambda, x1f, x2f)
end

h = 0.1 # time discretisation
tend = 200 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
xs[:,1] = [0.8; 0.83]
ys = zeros(N) # only one measurement

# Specify linear model
J = readcsv("j_ss1.csv")
ss = [0.856030685166100; 0.885965014362748] # steady state corresponding to J
lindu = begin
  A = eye(2)+h*J
  B = zeros(2,1)
  B[2] = h*delta
  b = -h*J*ss
  C = zeros(1,2)
  C[1] = 1.0
  D = zeros(1,1)
  Q = eye(2)*1e-7 # plant mismatch/noise
  R = eye(1)*1e-5 # measurement noise
  LLDS_functions.LLDS(A, B, b, C, D, Q, R)
end

linxs = zeros(2, N)
linys = zeros(N)
linxs[:, 1] = [0.8; 0.83]


# Simulate plant
ys[1] = (lindu.C*xs[:, 1] + rand(Normal(0.0, sqrt(lindu.R[1]))))[1] # measure from actual plant
linys[1] = (lindu.C*xs[:, 1] + rand(Normal(0.0, sqrt(lindu.R[1]))))[1] # model observations
for t=2:N
  xs[:, t] = DuReactor_functions.run_reactor(xs[:, t-1], 0.0, h, du) # actual plant
  ys[t] = (lindu.C*xs[:, t] + rand(Normal(0.0, sqrt(lindu.R[1]))))[1] # measured from actual plant
  xtemp, ytemp, ytemp_noiseless = LLDS_functions.step(linxs[:, t-1], [0.0], lindu)
  linys[t] = ytemp[1]
  linxs[:, t] = xtemp
end

figure(1) # Sanity check - the model and the plant coincide
x1, = plot(ts, xs[1,:]', "r")
x2, = plot(ts, xs[2,:]', "b")
y1, = plot(ts, ys, "r.")
linx1, = plot(ts, linxs[1,:]', "r--")
linx2, = plot(ts, linxs[2,:]', "b--")
liny1, = plot(ts, linys, "rx")
legend([x1,x2,y1,linx1,linx2, liny1],["Plant C_A","Plant T","Plant Measure C_A","Model C_A","Model T","Model Measure C_A"], loc="best")
