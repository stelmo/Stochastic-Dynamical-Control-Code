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
  C[2] = 0.0 # change this and get better results!
  D = zeros(1,1)
  Q = eye(2)*1e-6 # plant mismatch/noise
  R = eye(1)*1e-6 # measurement noise
  LLDS_functions.LLDS(A, B, b, C, D, Q, R)
end

linxs = zeros(2, N)
linys = zeros(N)
linxs[:, 1] = [0.8; 0.83]

us = zeros(N) # simulate some control movement. NOTE: us[1] = u(0), us[2] = u(1)... 
us[700:end] = 0.1
us[1200:end] = -0.1
# Simulate plant
ys[1] = (lindu.C*xs[:, 1] + rand(Normal(0.0, sqrt(lindu.R[1]))))[1] # measure from actual plant
linys[1] = (lindu.C*xs[:, 1] + rand(Normal(0.0, sqrt(lindu.R[1]))))[1] # model observations
for t=2:N
  xs[:, t] = DuReactor_functions.run_reactor(xs[:, t-1], us[t], h, du) # actual plant
  ys[t] = (lindu.C*xs[:, t] + rand(Normal(0.0, sqrt(lindu.R[1]))))[1] # measured from actual plant
  xtemp, ytemp, ytemp_noiseless = LLDS_functions.step(linxs[:, t-1], [us[t]], lindu)
  linys[t] = ytemp[1]
  linxs[:, t] = xtemp
end

# Filter
init_mean = [0.8; 0.83]
init_covar = eye(2)*1e-3 # vague
filtermeans = zeros(2, N)
filtercovars = zeros(2,2, N)
filtermeans[:, 1], filtercovars[:,:, 1] = LLDS_functions.init_filter(init_mean, init_covar, [us[1]], [ys[1]] , lindu)
for t=2:N
  filtermeans[:, t], filtercovars[:,:, t] = LLDS_functions.step_filter(filtermeans[:, t-1], filtercovars[:,:, t-1], [us[t]], [ys[t]], lindu)
end

# Prediction
pstart = 1000 # start predicting here (inclusive)
pend = N # prediction horizon
pred_us = zeros(1, pend-pstart+1)
pred_us[1,:] = us[pstart:pend]
pmeans, pcovars = LLDS_functions.predict_hidden(filtermeans[:, pstart-1], filtercovars[:,:, pstart-1], pred_us, lindu)


figure(1) # Sanity check - the model and the plant coincide (remember to set the control to 0)
suptitle("Modelling")
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "r")
linx1, = plot(ts, linxs[1,:]', "g--")
y1, = plot(ts, ys, "r.")
# liny1, = plot(ts, linys, "rx") # this is just to check that the noise order of magnitude is the same
legend([x1,y1,linx1],["Plant C_A","Plant Measure C_A","Model C_A"], loc="best")

subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "b")
linx2, = plot(ts, linxs[2,:]', "g--")
legend([x2,linx2],["Plant T","Model T"], loc="best")

figure(2) # check the filter results
suptitle("Filtering")
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "r")
k1, = plot(ts[1:10:end], filtermeans[1, 1:10:end]', "gx")
legend([x1,k1],["Plant C_A","Filter Mean C_A"], loc="best")

subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "b")
k2, = plot(ts[1:10:end], filtermeans[2, 1:10:end]', "gx")
legend([x2,k2],["Plant C_A","Filter Mean T"], loc="best")

figure(3) # check the prediction results
suptitle("Predicting")
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "r")
k1, = plot(ts[pstart:end], pmeans[1,:]', "gx")
legend([x1,k1],["Plant C_A","Predicted Mean C_A"], loc="best")

subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "b")
k2, = plot(ts[pstart:end], pmeans[2, :]', "gx")
legend([x2,k2],["Plant C_A","Predicted Mean T"], loc="best")
