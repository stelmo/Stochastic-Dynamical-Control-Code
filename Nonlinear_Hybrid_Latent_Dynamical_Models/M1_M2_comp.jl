# Implement the augmented switching dynamical system
using PyPlot
using Distributions
cd("..\\Linear_Hybrid_Latent_Dynamical_Models")
import SPF
reload("SPF.jl")
cd("..\\CSTR_Model")
using Reactor_functions
cd("..\\Linear_Latent_Dynamical_Models")
using Confidence
using LLDS_functions
cd("..\\Nonlinear_Hybrid_Latent_Dynamical_Models")

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
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
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
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

initial_states = [0.5; 450] # initial state
h = 0.1 # time discretisation
tend = 100.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
xsnofix = zeros(2, N)
ys1 = zeros(N)
ys2 = zeros(2, N)


C1 = [0.0 1.0]
R1 = eye(1)*10.0
C2 = eye(2)
R2 = eye(2)
R2[1] = 1e-3
R2[4] = 10.

Q = eye(2)
Q[1] = 1e-5
Q[4] = 4.0

A = [0.9 0.1;0.1 0.9]
# A = [0.5 0.5;0.5 0.5]

fun1(x,u,w) = Reactor_functions.run_reactor(x, u, h, cstr1)
fun2(x,u,w) = Reactor_functions.run_reactor(x, u, h, cstr2)
F = [fun1, fun2]

gs1(x) = C1*x
gs2(x) = C2*x
G1 = [gs1, gs1]
G2 = [gs2, gs2]

xdists = [MvNormal(Q); MvNormal(Q)]

ydists1 = [MvNormal(R1);MvNormal(R1)]
cstr_filter1 = SPF.Model(F, G1, A, xdists, ydists1)

ydists2 = [MvNormal(R2);MvNormal(R2)]
cstr_filter2 = SPF.Model(F, G2, A, xdists, ydists2)


nP = 500
initial_covar = eye(2)
initial_covar[1] = 1e-3
initial_covar[4] = 4.0
xdist = MvNormal(initial_states, initial_covar)
sdist = Categorical([0.5, 0.5])

particles1 = SPF.init_SPF(xdist, sdist, nP, 2)
fmeans1 = zeros(2, N)
fcovars1 = zeros(2, 2, N)

particles2 = SPF.init_SPF(xdist, sdist, nP, 2)
fmeans2 = zeros(2, N)
fcovars2 = zeros(2, 2, N)


us = zeros(N)

measurements1 = MvNormal(R1)
measurements2 = MvNormal(R2)
xs[:,1] = initial_states
xsnofix[:,1] = initial_states

ys1[1] = C1*xs[:, 1] + rand(measurements1) # measured from actual plant
SPF.init_filter!(particles1, 0.0, ys1[1], cstr_filter1)
fmeans1[:,1], fcovars1[:,:,1] = SPF.getStats(particles1)

ys2[:, 1] = C2*xs[:, 1] + rand(measurements2) # measured from actual plant
SPF.init_filter!(particles2, 0.0, ys2[:, 1], cstr_filter2)
fmeans2[:,1], fcovars2[:,:,1] = SPF.getStats(particles2)

# Loop through the rest of time
for t=2:N
  if ts[t] < 40.0
    xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], us[t-1], h, cstr1) # actual plant
    xsnofix[:, t] = Reactor_functions.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr1) # actual plant
  else
    xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], us[t-1], h, cstr2)
    xsnofix[:, t] = Reactor_functions.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr1) # actual plant
  end

  ys1[t] = C1*xs[:, t] + rand(measurements1) # measured from actual plant
  SPF.filter!(particles1, us[t-1], ys1[t], cstr_filter1)
  fmeans1[:,t], fcovars1[:,:,t] = SPF.getStats(particles1)

  ys2[:, t] = C2*xs[:, t] + rand(measurements2) # measured from actual plant
  SPF.filter!(particles2, us[t-1], ys2[:, t], cstr_filter2)
  fmeans2[:,t], fcovars2[:,:,t] = SPF.getStats(particles2)

end

rc("font", family="serif", size=24)

skip = 50
figure(2)
x1, = plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
x11, = plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
f1, = plot(fmeans1[1, 1:skip:end][:], fmeans1[2, 1:skip:end][:], "rx", markersize=5, markeredgewidth = 2)
f2, = plot(fmeans2[1, 1:skip:end][:], fmeans2[2, 1:skip:end][:], "bx", markersize=5, markeredgewidth = 2)
b1 = 0.0
b2 = 0.0
for k=1:skip:N
  p1, p2 = Confidence.plot95(fmeans1[:,k], fcovars1[:,:, k])
  b1, = plot(p1, p2, "r")

  p1, p2 = Confidence.plot95(fmeans2[:,k], fcovars2[:,:, k])
  b2, = plot(p1, p2, "b")
end
plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
ylabel("Temperature [K]")
xlabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1,f1, f2, b1, b2],["Nonlinear Model","Filter Mean: M1","Filter Mean: M2", L"Particle Filter M1: $1\sigma$-Ellipse",L"Particle Filter M2: $1\sigma$-Ellipse"], loc="best")
