# Implement the augmented switching dynamical system

using SPF
using Reactor
using Ellipse
using LLDS

# Add a definition for convert to make our lives easier!
# But be careful now!
function Base.convert(::Type{Float64}, x::Array{Float64, 1})
  return x[1]
end

# Specify the nonlinear model
cstr1 = begin
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

cstr2 = begin # slower reaction rate
  V = 5.0 #m3
  R = 8.314 #kJ/kmol.K
  CA0 = 1.0 #kmol/m3
  TA0 = 310.0 #K
  dH = -4.78e4 #kJ/kmol
  k0 = 72.0e6 #1/min # changed here!
  E = 8.314e4 #kJ/kmol
  Cp = 0.239 #kJ/kgK # changed here!
  rho = 1000.0 #kg/m3
  F = 100e-3 #m3/min
  Reactor.reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

initial_states = [0.5; 450] # initial state
h = 0.1 # time discretisation
tend = 150.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
xsnofix = zeros(2, N)
ys = zeros(2, N)

newC = eye(2)
R = eye(2)
R[1] = 1e-3
R[4] = 10.
Q = eye(2)
Q[1] = 1e-6
Q[4] = 0.01

A = [0.9 0.1;0.1 0.9]
# A = [0.5 0.5;0.5 0.5]
fun1(x,u,w) = Reactor.run_reactor(x, u, h, cstr1)
fun2(x,u,w) = Reactor.run_reactor(x, u, h, cstr2)
gs(x) = newC*x
F = [fun1, fun2]
G = [gs, gs]

ydists = [MvNormal(R);MvNormal(R)]
xdists = [MvNormal(Q); MvNormal(Q)]
cstr_filter = SPF.Model(F, G, A, xdists, ydists)

nP = 1500
initial_covar = eye(2)
initial_covar[1] = 1e-3
initial_covar[4] = 4.0
xdist = MvNormal(initial_states, initial_covar)
sdist = Categorical([0.5, 0.5])
particles = SPF.init_SPF(xdist, sdist, nP, 2)
fmeans = zeros(2, N)
fcovars = zeros(2,2, N)

us = zeros(N)
switchtrack = zeros(2, N)
state_dist = MvNormal(Q)

measurements = MvNormal(R)
xs[:,1] = initial_states
xsnofix[:,1] = initial_states

ys[:, 1] = newC*xs[:, 1] + rand(measurements) # measured from actual plant
SPF.init_filter!(particles, 0.0, ys[:, 1], cstr_filter)

for k=1:2
  switchtrack[k, 1] = sum(particles.w[find((x)->x==k, particles.s)])
end

fmeans[:,1], fcovars[:,:,1] = SPF.getStats(particles)
# Loop through the rest of time
for t=2:N
  random_element = rand(state_dist)
  if ts[t] < 50.0
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr1) + random_element # actual plant
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr1) + random_element # actual plant
  else
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr2) + random_element
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr1) + random_element # actual plant
  end

  ys[:, t] = newC*xs[:, t] + rand(measurements) # measured from actual plant
  SPF.filter!(particles, us[t-1], ys[:, t], cstr_filter)
  fmeans[:,t], fcovars[:,:,t] = SPF.getStats(particles)

  for k=1:2
    switchtrack[k, t] = sum(particles.w[find((x)->x==k, particles.s)])
  end

end

rc("font", family="serif", size=24)

figure(1)
maxswitch = maximum(switchtrack)
axes = Array(Any, 2)
im = 0
width = 500
for k=1:2
  ax = subplot(2, 1, k)
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


skip = 50
figure(2)
x1, = plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
x11, = plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
f1, = plot(fmeans[1, 1:skip:end][:], fmeans[2, 1:skip:end][:], "rx", markersize=5, markeredgewidth = 2)
b1 = 0.0
for k=1:skip:N
  p1, p2 = Ellipse.ellipse(fmeans[:,k], fcovars[:,:, k])
  b1, = plot(p1, p2, "b")
end
plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
ylabel("Temperature [K]")
xlabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1,f1, b1],["Nonlinear Model","Particle Filter Mean", L"Particle Filter $1\sigma$-Ellipse"], loc="best")

skipm = skip
figure(3) # Plot filtered results
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
x1nf, = plot(ts, xsnofix[1,:]', "g--", linewidth=3)
y2, = plot(ts[1:skipm:end], ys[1, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k1, = plot(ts, fmeans[1,:]', "r--", linewidth=3)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1, k1],["Nonlinear Model","Filtered Mean"], loc="best")
xlim([0, tend])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
x2nf, = plot(ts, xsnofix[2,:]', "g--", linewidth=3)
y2, = plot(ts[1:skipm:end], ys[2, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts, fmeans[2,:]', "r--", linewidth=3)
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y2, x2nf],["Nonlinear Model Measured","Nonlinear Model No Switch"], loc="best")
xlim([0, tend])
