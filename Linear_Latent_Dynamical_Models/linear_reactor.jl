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

init_state = [0.57; 395]
h = 0.001 # time discretisation
tend = 1.0 # end simulation time
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
  Q[1] = 2.5e-7
  Q[4] = 1e-3
  R = 2.0 # measurement noise
  LLDS_functions.LLDS{Float64}(A, B, b, C, Q, R)
end

# Plant initialisation
linxs = zeros(2, N)
linys = zeros(N)
linxs[:, 1] = init_state
us = zeros(N) # simulate some control movement. NOTE: us[1] = u(t=0), us[2] =u(t=1)...

# Filter
init_mean = init_state
init_covar = eye(2) # vague
init_covar[1] = 1e-6
init_covar[4] = 1e-6
filtermeans = zeros(2, N)
filtercovars = zeros(2,2, N)
filtermeans[:, 1], filtercovars[:,:, 1] = LLDS_functions.init_filter(init_mean, init_covar, ys[1], lin_cstr)

# Simulate plant
norm_dist = Normal(0.0, sqrt(lin_cstr.R))
ys[1] = lin_cstr.C*xs[:, 1] + rand(norm_dist) # measure from actual plant
for t=2:N
  xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], us[t], h, cstr) # actual plant
  ys[t] = lin_cstr.C*xs[:, t] + rand(norm_dist) # measured from actual plant
  linxs[:, t], temp1, temp2 = LLDS_functions.step(linxs[:, t-1], us[t], lin_cstr)
  filtermeans[:, t], filtercovars[:,:, t] = LLDS_functions.step_filter(filtermeans[:, t-1], filtercovars[:,:, t-1], us[t], ys[t], lin_cstr)
end

# Prediction
# pstart = 100 # start predicting here (inclusive)
# pend = N # prediction horizon
# pred_us = zeros(pend-pstart+1)
# pred_us[:] = us[pstart-1:pend-1]
# pmeans, pcovars = LLDS_functions.predict_hidden(filtermeans[:, pstart-1], filtercovars[:,:, pstart-1], pred_us, lin_cstr)


# figure(1) # Sanity check - the model and the plant coincide (remember to set the control to 0)
# subplot(2,1,1)
# x1, = plot(ts, xs[1,:]', "k", linewidth=3)
# linx1, = plot(ts, linxs[1,:]', "r--", linewidth=3)
# ylabel(L"Concentration $[kmol.m^{-3}]$")
# legend([x1,linx1],[L"Nonlinear Model $C_A$",L"Linear Model $C_A$"], loc="best")
#
# subplot(2,1,2)
# x2, = plot(ts, xs[2,:]', "k", linewidth=3)
# linx2, = plot(ts, linxs[2,:]', "r--", linewidth=3)
# y2, = plot(ts, ys, "rx", markersize=4)
# ylabel(L"Temperature $[K]$")
# xlabel(L"Time $[min]$")
# legend([x2, linx2, y2],[L"Nonlinear Model $T_R$",L"Linear Model $T_R$", L"Nonlinear Model Measured $T_R$"], loc="best")
#
# rc("font",size=22)
#
# figure(2) # check the filter results
# subplot(2,1,1)
# x1, = plot(ts, xs[1,:]', "k")
# k1, = plot(ts[1:10:end], filtermeans[1, 1:10:end]', "mo")
# ylabel(L"Concentration $[kmol.m^{-3}]$")
# legend([x1,k1],[L"Nonlinear Model $C_A$",L"Filtered $C_A$"], loc="best")
#
# subplot(2,1,2)
# x2, = plot(ts, xs[2,:]', "k")
# k2, = plot(ts[1:10:end], filtermeans[2, 1:10:end]', "mo")
# y2, = plot(ts, ys, "rx", markersize=4)
# ylabel(L"Temperature $[K]$")
# xlabel(L"Time $[min]$")
# legend([x2,k2, y2],[L"Nonlinear Model $T_R$",L"Filtered $T_R$",L"Nonlinear Model Measured $T_R$"], loc="best")
#
# rc("font",size=22)

figure(3) # Combined
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
linx1, = plot(ts, linxs[1,:]', "r--", linewidth=3)
k1, = plot(ts[1:10:end], filtermeans[1, 1:10:end]', "mo")
ylabel(L"Concentration $[kmol.m^{-3}]$")
legend([x1,linx1, k1],["Nonlinear Model","Noisy Linear Model", "Filtered"], loc="best")
xlim([0, tend])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
linx2, = plot(ts, linxs[2,:]', "r--", linewidth=3)
y2, = plot(ts, ys, "rx", markersize=5, markeredgewidth=1)
k2, = plot(ts[1:10:end], filtermeans[2, 1:10:end]', "mo")
ylabel(L"Temperature $[K]$")
xlabel(L"Time $[min]$")
legend([y2],["Nonlinear Model Measured"], loc="best")
xlim([0, tend])
rc("font",size=22)



# figure(3) # check the prediction results
# suptitle("Predicting")
# subplot(2,1,1)
# x1, = plot(ts, xs[1,:]', "r")
# k1, = plot(ts[pstart:end], pmeans[1,:]', "gx")
# legend([x1,k1],[L"Nonlinear Model $C_A$",L"Predicted $C_A$"], loc="best")
#
# subplot(2,1,2)
# x2, = plot(ts, xs[2,:]', "b")
# k2, = plot(ts[pstart:end], pmeans[2, :]', "gx")
# legend([x2,k2],[L"Nonlinear Model $T_R$",L"Predicted $T_R$"], loc="best")
