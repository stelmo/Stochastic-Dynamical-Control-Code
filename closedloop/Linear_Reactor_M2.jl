# This script will run the linear cstr reactor.

using LLDS
using Reactor
using Ellipse
using LQR

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

init_state = [0.5; 450]
h = 0.1 # time discretisation
tend = 2000.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
xs[:,1] = init_state
ys = zeros(2, N) # only one measurement

xspace = [0.0, 1.0]
yspace = [250, 550]

# Specify the linear model
linsystems = Reactor.getLinearSystems_randomly(0, xspace, yspace, h, cstr) # doesnt work weirdly...
opoint = 2
A = linsystems[opoint].A
B = linsystems[opoint].B
b = linsystems[opoint].b
C = eye(2)
Q = eye(2) # plant mismatch/noise
Q[1] = 1e-10 #1e-06
Q[4] = 0.00001 #0.1
R = eye(2)
R[1] = 1e-3
R[4] = 10.0 # measurement noise
lin_cstr = LLDS.llds(A, B, C, Q, R) # so that the KF works

# Controller
QQ = zeros(2, 2)
QQ[1] = 1.0
RR = 1.0
H = [1.0 0.0]
ysp = linsystems[opoint].op[1] - b[1] # remember to adjust for ss standard
# ysp = 0.01 - b[1]
# ysp = linsystems[3].op[1] - b[1] # remember to adjust for ss standard
x_off, u_off = LQR.offset(A,B,C,H, ysp)
K = LQR.lqr(A, B, QQ, RR)

# Plant initialisation
us = zeros(N-1) # simulate some control movement.

# Simulate plant
state_dist = MvNormal(Q)
norm_dist = MvNormal(lin_cstr.R)
ys[:, 1] = C*xs[:, 1] + rand(norm_dist) # measure from actual plant
# Filter
init_mean = init_state - b
init_covar = eye(2) # vague
init_covar[1] = 1e-3
init_covar[4] = 4.
filtermeans = zeros(2, N)
filtercovars = zeros(2,2, N)
filtermeans[:, 1], filtercovars[:,:, 1] = LLDS.init_filter(init_mean, init_covar, ys[:, 1]-b, lin_cstr)

for t=2:N
  if t%10 == 0 || t==2 # to start
    us[t-1] = -K*(filtermeans[:, t-1] - x_off) + u_off # controller action
  else
    us[t-1] = us[t-2]
  end
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr) + rand(state_dist) # actual plant
  ys[:, t] = lin_cstr.C*xs[:, t] + rand(norm_dist) # measured from actual plant
  filtermeans[:, t], filtercovars[:,:, t] = LLDS.step_filter(filtermeans[:, t-1], filtercovars[:,:, t-1], us[t-1], ys[:, t]-b, lin_cstr)
end
filtermeans = filtermeans .+ b


rc("font", family="serif", size=24)


skipmeas = int(length(ts)/20)
skipmean = int(length(ts)/20)
figure(1) # Filtering
subplot(3,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
y2, = plot(ts[1:skipmeas:end], ys[1, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
k1, = plot(ts[1:skipmean:end], filtermeans[1, 1:skipmean:end]', "bx", markersize=5, markeredgewidth = 2)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([y2, x1],["Nonlinear Model Measured", "Nonlinear Model"], loc="best")
xlim([0, tend])
ylim([0, 1])

subplot(3,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
y2, = plot(ts[1:skipmeas:end], ys[2, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts[1:skipmean:end], filtermeans[2, 1:skipmean:end]', "bx", markersize=5, markeredgewidth = 2)
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([k2],["Filtered Mean Estimate"], loc="best")
xlim([0, tend])
ylim([minimum(xs[2,:]), maximum(xs[2,:])])
subplot(3,1,3)
plot(ts[1:end-1], us)
