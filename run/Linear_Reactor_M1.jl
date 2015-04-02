# This script will run the linear cstr reactor.

using LLDS
using Reactor
using Ellipse

# This is to add the noise to the measurements! cant add a float to a [float]...
function Base.convert(::Type{Float64}, x::Array{Float64, 1})
  return x[1]
end

# Specify the nonlinear model
cstr = begin
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

init_state = [0.50; 400]
h = 0.1 # time discretisation
tend = 150.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
xs[:,1] = init_state
ys = zeros(N) # only one measurement

xspace = [0.0, 1.0]
yspace = [250, 550]

# Specify the linear model
npoints = 1
linsystems = Reactor.getLinearSystems_randomly(npoints, xspace, yspace, h, cstr) # doesnt work weirdly...
A = linsystems[3].A
B = linsystems[3].B
b = linsystems[3].b
C = [0.0 1.0]
Q = eye(2) # plant mismatch/noise
Q[1] = 1e-5
Q[4] = 4.
R = 10.0 # measurement noise
lin_cstr = LLDS.llds(A, B, b, C, Q, R)


# Plant initialisation
linxs = zeros(2, N)
linxs[:, 1] = init_state
us = zeros(N) # simulate some control movement. NOTE: us[1] = u(t=0), us[2] =u(t=1)...


# Simulate plant
state_dist = MvNormal(Q)
norm_dist = Normal(0.0, sqrt(lin_cstr.R))
ys[1] = lin_cstr.C*xs[:, 1] + rand(norm_dist) # measure from actual plant
# Filter
init_mean = init_state
init_covar = eye(2) # vague
init_covar[1] = 1e-3
init_covar[4] = 10.
filtermeans = zeros(2, N)
filtercovars = zeros(2,2, N)
filtermeans[:, 1], filtercovars[:,:, 1] = LLDS.init_filter(init_mean, init_covar, ys[1], lin_cstr)
for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t], h, cstr) + rand(state_dist)# actual plant
  ys[t] = lin_cstr.C*xs[:, t] + rand(norm_dist) # measured from actual plant
  linxs[:, t], temp = LLDS.step(linxs[:, t-1], us[t], lin_cstr)
  filtermeans[:, t], filtercovars[:,:, t] = LLDS.step_filter(filtermeans[:, t-1], filtercovars[:,:, t-1], us[t], ys[t], lin_cstr)
end

# Prediction
# pstart = 100 # start predicting here (inclusive)
# pend = N # prediction horizon
# pred_us = zeros(pend-pstart+1)
# pred_us[:] = us[pstart-1:pend-1]
# pmeans, pcovars = LLDS_functions.predict_hidden(filtermeans[:, pstart-1], filtercovars[:,:, pstart-1], pred_us, lin_cstr)

rc("font", family="serif", size=24)

skip = 120
figure(1)
x1, = plot(xs[1,:][:], xs[2,:][:], "k",linewidth=3)
f1, = plot(filtermeans[1, 1:skip:end][:], filtermeans[2, 1:skip:end][:], "rx", markersize=5, markeredgewidth = 2)
b1 = 0.0
for k=1:skip:N
  p1, p2 = Ellipse.ellipse(filtermeans[:,k], filtercovars[:,:, k])
  b1, = plot(p1, p2, "b")
end
plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
ylabel("Temperature [K]")
xlabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1,f1, b1],["Nonlinear Model","Kalman Filter Mean", L"Kalman Filter $1\sigma$-Ellipse"], loc="best")


skipm = skip
figure(2) # Filtering
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
linx1, = plot(ts, linxs[1,:]', "r--", linewidth=3)
k1, = plot(ts[1:skip:end], filtermeans[1, 1:skip:end]', "mo")
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1,linx1],["Nonlinear Model","Linear Model"], loc="best")
xlim([0, tend])
ylim([0, 1.5])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
linx2, = plot(ts, linxs[2,:]', "r--", linewidth=3)
y2, = plot(ts[1:skipm:end], ys[1:skipm:end], "rx", markersize=5, markeredgewidth=1)
k2, = plot(ts[1:skip:end], filtermeans[2, 1:skip:end]', "mo")
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y2, k2],["Nonlinear Model Measured", "Filtered Mean Estimate"], loc="best")
xlim([0, tend])
ylim([minimum(xs[2,:]), maximum(xs[2,:])])
