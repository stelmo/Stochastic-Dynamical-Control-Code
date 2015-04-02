# Implement the augmented switching dynamical system
using PyPlot
using Distributions

import RBPF
import SPF
cd("..\\CSTR_Model")
using Reactor_functions
cd("..\\Linear_Latent_Dynamical_Models")
using Confidence
using LLDS_functions
cd("..\\Linear_Hybrid_Latent_Dynamical_Models")

# Add a definition for convert to make our lives easier!
# But be careful now!
function Base.convert(::Type{Float64}, x::Array{Float64, 1})
  return x[1]
end

# Specify the nonlinear model
cstr_model = begin
  V = 5.0; #m3
  R = 8.314; #kJ/kmol.K
  CA0 = 1.0; #kmol/m3
  TA0 = 310.0; #K
  dH = -4.78e4; #kJ/kmol
  k0 = 72.0e7; #1/min
  E = 8.314e4; #kJ/kmol
  Cp = 0.239; #kJ/kgK
  rho = 1000.0; #kg/m3
  F = 100e-3; #m3/min
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

h = 0.1 # time discretisation
tend = 150. # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
ys = zeros(2, N) # only one measurement

init_state = [0.5; 400] # initial state
C = eye(2) # observe both states
R = eye(2)
R[1] = 1e-3
R[4] = 10.0
Q = eye(2)
Q[1] = 1e-5
Q[4] = 4.0

# Divide state space into sectors: n by m
nX = 2 # rows
nY = 2 # cols
xspace = [0.0, 1.0]
yspace = [250, 650]

linsystems = Reactor_functions.getLinearSystems_randomly(0, xspace, yspace, h, cstr_model)
A = linsystems[2].A
B = linsystems[2].B
b = linsystems[2].b
lin_cstr = LLDS_functions.LLDS(A, B, b, C, Q, R)


linsystems = Reactor_functions.getLinearSystems(nX, nY, xspace, yspace, h, cstr_model)
# linsystems = Reactor_functions.getLinearSystems_randomly(0, xspace, yspace, h, cstr_model)

models, A = RBPF.setup_RBPF(linsystems, C, Q, R)

nP = 500
initial_states = init_state
initial_covar = eye(2)
initial_covar[1] = 1e-3
initial_covar[4] = 4.0
sguess =  SPF.getInitialSwitches(initial_states, linsystems)
particles = RBPF.init_RBPF(Categorical(sguess), initial_states, initial_covar, 2, nP)

fmeans = zeros(2, N)
fcovars = zeros(2,2,N)
filtermeans = zeros(2, N)
filtercovars = zeros(2,2, N)
switchtrack = zeros(length(linsystems), N)

state_dist = MvNormal(Q)
us = zeros(N)
measurements = MvNormal(R)
xs[:,1] = initial_states
ys[:, 1] = C*xs[:, 1] + rand(measurements) # measured from actual plant
RBPF.init_filter!(particles, 0.0, ys[:, 1], models)
fmeans[:,1], fcovars[:,:, 1] = RBPF.getStats(particles)

filtermeans[:, 1], filtercovars[:,:, 1] = LLDS_functions.init_filter(initial_states, initial_covar, ys[:, 1], lin_cstr)

for k=1:length(linsystems)
  switchtrack[k, 1] = sum(particles.ws[find((x)->x==k, particles.ss)])
end
# Loop through the rest of time
for t=2:N
  xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_dist) # actual plant
  ys[:, t] = C*xs[:, t] + rand(measurements) # measured from actual plant
  RBPF.filter!(particles, us[t-1], ys[:, t], models, A)
  fmeans[:, t], fcovars[:,:, t] = RBPF.getStats(particles)
  for k=1:length(linsystems)
    switchtrack[k, t] = sum(particles.ws[find((x)->x==k, particles.ss)])
  end
  filtermeans[:, t], filtercovars[:,:, t] = LLDS_functions.step_filter(filtermeans[:, t-1], filtercovars[:,:, t-1], us[t], ys[:, t], lin_cstr)
end

rc("font", family="serif", size=24)

figure(1)
for k=1:length(linsystems)
  plot(linsystems[k].op[1],linsystems[k].op[2],"kx",markersize=5, markeredgewidth=1)
annotate(string("Switch: ", k),
      xy=[linsystems[k].op[1],linsystems[k].op[2]],
      xytext=[linsystems[k].op[1],linsystems[k].op[2]],
      fontsize=22.0,
      ha="center",
      va="bottom")
end
plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
xlim([-0.1, 1.1])
xlabel(L"Concentration [kmol.m$^{-3}$]")
ylabel("Temperature [K]")

figure(2)
maxswitch = maximum(switchtrack)
axes = Array(Any, length(linsystems))
im = 0
width = 500
for k=1:length(linsystems)
  ax = subplot(length(linsystems), 1, k)
  axes[k] = ax
  im = imshow(repeat(switchtrack[k,:], outer=[width, 1]), cmap="cubehelix",vmin=0.0, vmax=maxswitch, interpolation="nearest", aspect="auto")
  tick_params(axis="y", which="both",left="off",right="off", labelleft = "off")
  tick_params(axis="x", which="both",bottom="off", labelbottom = "off")
  ylabel(string("S::",k))
end
tick_params(axis="x", labelbottom = "on")
xticks([1:int(length(ts)/10.0):length(ts)], ts[1:int(length(ts)/10.0):end])
colorbar(im, ax=axes)
xlabel("Time [min]")

figure(3) # Plot filtered results
skip = int(length(ts)/75)
skipm = skip
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
y2, = plot(ts[1:skipm:end], ys[1, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k1, = plot(ts[1:skip:end], fmeans[1, 1:skip:end]', "r--",linewidth=3)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1],["Nonlinear Model"], loc="best")
xlim([0, tend])
ylim([0, 1])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
y2, = plot(ts[1:skipm:end], ys[2, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts[1:skip:end], fmeans[2, 1:skip:end]', "r--",linewidth=3)
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y2, k2],["Nonlinear Model Measured", "Filtered Mean Estimate"], loc="best")
xlim([0, tend])


skip = int(length(ts)/20)
figure(4) # Kalman Filter Demonstration
x1, = plot(xs[1,:][:], xs[2,:][:], "k",linewidth=3)
x11, = plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
f1, = plot(fmeans[1, 1:skip:end][:], fmeans[2, 1:skip:end][:], "rx", markersize=5, markeredgewidth = 2)
f2, = plot(filtermeans[1, 1:skip:end][:], filtermeans[2, 1:skip:end][:], "gx", markersize=5, markeredgewidth = 2)
b1 = 0.0
b2 = 0.0
for k=1:skip:N
  p1, p2 = Confidence.plot95(fmeans[:,k], fcovars[:,:, k])
  b1, = plot(p1, p2, "b")

  p3, p4 = Confidence.plot95(filtermeans[:,k], filtercovars[:,:, k])
  b2, = plot(p3, p4, "g")

end
plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
ylabel("Temperature [K]")
xlabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1,f1,f2, b1, b2],["Nonlinear Model","Switching Kalman Filter Mean","Kalman Filter Mean", L"Switching Kalman Filter $1\sigma$-Ellipse",L"Kalman Filter $1\sigma$-Ellipse"], loc="best")
