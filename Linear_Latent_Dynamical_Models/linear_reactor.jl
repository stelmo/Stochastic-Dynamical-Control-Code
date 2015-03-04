# This script will run the linear cstr reactor.

using PyPlot
using Distributions
import LLDS_functions
import Reactor_functions

# Add a definition for convert to make our lives easier!
# But be careful now!
function Base.convert(::Type{Float64}, x::Array{Float64, 1})
  return x[1]
end
function Base.convert(::Type{Float64}, x::Array{Float64, 2})
  return x[1]
end

# Specify the nonlinear model
cstr = begin
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

init_state = [0.05; 520]
h = 0.001 # time discretisation
tend = 0.5 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
xs[:,1] = init_state
ys = zeros(N) # only one measurement

# Specify the linear model
J = readcsv("J_ss.csv")
ss = readcsv("ss.csv")'[:,1]
lin_cstr = begin
  A = eye(2)+h*J
  B = zeros(2,1)
  B[2] = 1./(cstr.rho*cstr.Cp*cstr.V)
  b = -h*J*ss
  C = zeros(1,2)
  C[2] = 1.0 #measure temperature
  Q = eye(2) # plant mismatch/noise
  Q[1] = 1e-6
  Q[4] = 1e-10
  R = 2.0 # measurement noise
  LLDS_functions.LLDS{Float64}(A, B, b, C, Q, R)
end

linxs = zeros(2, N)
linys = zeros(N)
linxs[:, 1] = init_state

us = zeros(N) # simulate some control movement. NOTE: us[1] = u(1), us[2] = u(2)...
# There is no control
# us[700:end] = 0.0
# us[1200:end] = 0.0
# Simulate plant
norm_dist = Normal(0.0, sqrt(lin_cstr.R))
ys[1] = lin_cstr.C*xs[:, 1] + rand(norm_dist) # measure from actual plant
linys[1] = lin_cstr.C*xs[:, 1] + rand(norm_dist) # model observations
for t=2:N
  xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], us[t], h, cstr) # actual plant
  ys[t] = lin_cstr.C*xs[:, t] + rand(norm_dist) # measured from actual plant
  linxs[:, t], linys[t], ytemp_noiseless = LLDS_functions.step(linxs[:, t-1], us[t], lin_cstr)
end
#
# # Filter
# init_mean = init_state
# init_covar = eye(2)*1e-3 # vague
# filtermeans = zeros(2, N)
# filtercovars = zeros(2,2, N)
# filtermeans[:, 1], filtercovars[:,:, 1] = LLDS_functions.init_filter(init_mean, init_covar, ys[1], lin_cstr)
# for t=2:N
#   filtermeans[:, t], filtercovars[:,:, t] = LLDS_functions.step_filter(filtermeans[:, t-1], filtercovars[:,:, t-1], us[t], ys[t], lin_cstr)
# end
#
# # Prediction
# pstart = 1000 # start predicting here (inclusive)
# pend = N # prediction horizon
# pred_us = zeros(pend-pstart+1)
# pred_us[:] = us[pstart-1:pend-1]
# pmeans, pcovars = LLDS_functions.predict_hidden(filtermeans[:, pstart-1], filtercovars[:,:, pstart-1], pred_us, lin_cstr)


figure(1) # Sanity check - the model and the plant coincide (remember to set the control to 0)
suptitle("Modelling")
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "r")
linx1, = plot(ts, linxs[1,:]', "g--")
# liny1, = plot(ts, linys, "rx") # this is just to check that the noise order of magnitude is the same
legend([x1,y1,linx1],[L"Nonlinear Model $C_A$",L"Linear Model $C_A$"], loc="best")

subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "b")
linx2, = plot(ts, linxs[2,:]', "g--")
y2, = plot(ts, ys, "b.")
legend([x2,linx2, y2],[L"Nonlinear Model $T_R$",L"Linear Model $T_R$", L"Nonlinear Model Measured $T_R$"], loc="best")
#
# figure(2) # check the filter results
# suptitle("Filtering")
# subplot(2,1,1)
# x1, = plot(ts, xs[1,:]', "r")
# k1, = plot(ts[1:10:end], filtermeans[1, 1:10:end]', "gx")
# legend([x1,k1],["Plant C_A","Filter Mean C_A"], loc="best")
#
# subplot(2,1,2)
# x2, = plot(ts, xs[2,:]', "b")
# k2, = plot(ts[1:10:end], filtermeans[2, 1:10:end]', "gx")
# legend([x2,k2],["Plant C_A","Filter Mean T"], loc="best")
#
# figure(3) # check the prediction results
# suptitle("Predicting")
# subplot(2,1,1)
# x1, = plot(ts, xs[1,:]', "r")
# k1, = plot(ts[pstart:end], pmeans[1,:]', "gx")
# legend([x1,k1],["Plant C_A","Predicted Mean C_A"], loc="best")
#
# subplot(2,1,2)
# x2, = plot(ts, xs[2,:]', "b")
# k2, = plot(ts[pstart:end], pmeans[2, :]', "gx")
# legend([x2,k2],["Plant C_A","Predicted Mean T"], loc="best")
