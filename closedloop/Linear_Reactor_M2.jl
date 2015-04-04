# This script will run the linear cstr reactor.

using LLDS
using Reactor
using Ellipse


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
tend = 20.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
xs[:,1] = init_state
ys = zeros(2, N)

xspace = [0.0, 1.0]
yspace = [250, 550]

linsystems = Reactor.getLinearSystems_randomly(0, xspace, yspace, h, cstr)
A = linsystems[2].A
B = linsystems[2].B
b = linsystems[2].b
C = eye(2)
Q = eye(2) # plant mismatch/noise
Q[1] = 1e-5
Q[4] = 4.
R = eye(2)
R[1] = 1e-3
R[4] = 10.0 # measurement noise
lin_cstr = LLDS.llds(A, B, b, C, Q, R)


# Plant initialisation
us = zeros(N) # simulate some control movement. NOTE: us[1] = u(t=0), us[2] =u(t=1)...


# Simulate plant
state_dist = MvNormal(Q)
norm_dist = MvNormal(lin_cstr.R)
ys[:, 1] = lin_cstr.C*xs[:, 1] + rand(norm_dist) # measure from actual plant
# Filter
init_mean = init_state
init_covar = eye(2) # vague
init_covar[1] = 1e-3
init_covar[4] = 4.
filtermeans = zeros(2, N)
filtercovars = zeros(2,2, N)
filtermeans[:, 1], filtercovars[:,:, 1] = LLDS.init_filter(init_mean, init_covar, ys[:, 1], lin_cstr)
for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t], h, cstr) + rand(state_dist)
  ys[:, t] = lin_cstr.C*xs[:, t] + rand(norm_dist) # measured from actual plant
  filtermeans[:, t], filtercovars[:,:, t] = LLDS.step_filter(filtermeans[:, t-1], filtercovars[:,:, t-1], us[t-1], ys[:, t], lin_cstr)
end

# Prediction
pstart = 100 # start predicting here (inclusive)
pend = N # prediction horizon
pred_us = zeros(pend-pstart+1)
pred_us = us[pstart-1:pend-1]
pmeans, pcovars = LLDS.predict_hidden(filtermeans[:, pstart-1], filtercovars[:,:, pstart-1], pred_us, lin_cstr)

rc("font", family="serif", size=24)

skip = 20
figure(1) # Kalman Filter Demonstration
x1, = plot(xs[1,:][:], xs[2,:][:], "k",linewidth=3)
f1, = plot(filtermeans[1, 1:skip:end][:], filtermeans[2, 1:skip:end][:], "bx", markersize=5, markeredgewidth = 2)
b1 = 0.0
for k=1:skip:N
  p1, p2 = Ellipse.ellipse(filtermeans[:,k], filtercovars[:,:, k])
  b1, = plot(p1, p2, "b")
end
np1, np2 = size(pmeans)
f1, = plot(pmeans[1, 1:skip:end][:], pmeans[2, 1:skip:end][:], "rx", markersize=5, markeredgewidth = 2)
b2 = 0.0
for k=1:skip:np2
  p1, p2 = Ellipse.ellipse(pmeans[:,k], pcovars[:,:, k])
  b2, = plot(p1, p2, "r")
end
plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
ylabel("Temperature [K]")
xlabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1,f1, b1, b2],["Nonlinear Model","Kalman Filter Mean", L"Kalman Filter $1\sigma$-Ellipse","Prediction"], loc="best")
