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
  k0 = 72.0e8 #1/min # changed here!
  E = 8.314e4 #kJ/kmol
  Cp = 0.239 #kJ/kgK # changed here!
  rho = 1000.0 #kg/m3
  F = 100e-3 #m3/min
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

initial_states = [0.5; 450] # initial state
h = 0.001 # time discretisation
tend = 5.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
xsnofix = zeros(2, N)
ys = zeros(N)

newC = [0.0 1.0]
R = eye(1)*4.0
Q = eye(2)
Q[1] = 1e-4
Q[4] = 5.0

# A = [0.8 0.2;0.2 0.8]
A = [0.5 0.5;0.5 0.5]
fun1(x,u,w) = Reactor_functions.run_reactor(x, u, h, cstr1)
fun2(x,u,w) = Reactor_functions.run_reactor(x, u, h, cstr2)
gs(x) = newC*x
F = [fun1, fun2]
G = [gs, gs]

ydists = [MvNormal(R);MvNormal(R)]
xdists = [MvNormal(Q); MvNormal(Q)]
cstr_filter = SPF.Model(F, G, A, xdists, ydists)

nP = 500
initial_covar = eye(2)
initial_covar[1] = 1e-6
initial_covar[4] = 1.0
xdist = MvNormal(initial_states, initial_covar)
sdist = Categorical([0.5, 0.5])
particles = SPF.init_SPF(xdist, sdist, nP, 2)
fmeans = zeros(2, N)
fcovars = zeros(2,2, N)

us = zeros(N)
switch_count = zeros(2, N)

measurements = MvNormal(R)
xs[:,1] = initial_states
xsnofix[:,1] = initial_states

ys[1] = newC*xs[:, 1] + rand(measurements) # measured from actual plant
SPF.init_filter!(particles, 0.0, ys[1], cstr_filter)

for k=1:nP
  if particles.s[k] == 1
    switch_count[1, 1] += particles.w[k]
  else
    switch_count[2, 1] += particles.w[k]
  end
end

fmeans[:,1], fcovars[:,:,1] = SPF.getStats(particles)
# Loop through the rest of time
for t=2:N
  if ts[t] < 1.0
    xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], us[t-1], h, cstr1) # actual plant
    xsnofix[:, t] = Reactor_functions.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr1) # actual plant
  else
    xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], us[t-1], h, cstr2)
    xsnofix[:, t] = Reactor_functions.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr1) # actual plant
  end
  # if ts[t] > 0.5
  #   us[t] = -300.0
  # end
  ys[t] = newC*xs[:, t] + rand(measurements) # measured from actual plant
  SPF.filter!(particles, us[t-1], ys[t], cstr_filter)
  fmeans[:,t], fcovars[:,:,t] = SPF.getStats(particles)

  for k=1:nP
    if particles.s[k] == 1
      switch_count[1, t] += particles.w[k]
    else
      switch_count[2, t] += particles.w[k]
    end
  end

end


figure(1)
plot(ts, switch_count[1,:][:], "g")
plot(ts, switch_count[2,:][:], "k")
xlabel("Time [min]")
ylabel("Weight")
rc("font",size=22)


skip = 50
# figure(2) # Kalman Filter Demonstration
# x1, = plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
# f1, = plot(fmeans[1, 1:skip:end][:], fmeans[2, 1:skip:end][:], "rx", markersize=5, markeredgewidth = 2)
# b1 = 0.0
# for k=1:skip:N
#   p1, p2 = Confidence.plot95(fmeans[:,k], fcovars[:,:, k])
#   b1, = plot(p1, p2, "b")
# end
# ylabel("Temperature [K]")
# xlabel(L"Concentration [kmol.m$^{-3}$]")
# legend([x1,f1, b1],["Nonlinear Model","Particle Filter Mean", L"Particle Filter $1\sigma$-Ellipse"], loc="best")


figure(3) # Plot filtered results
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
x1nf, = plot(ts, xsnofix[1,:]', "g--", linewidth=3)
k1, = plot(ts, fmeans[1,:]', "r--", linewidth=3)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1, k1],["Nonlinear Model","Filtered Mean"], loc="best")
xlim([0, tend])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
x2nf, = plot(ts, xsnofix[2,:]', "g--", linewidth=3)
y2, = plot(ts[1:10:end], ys[1:10:end], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts, fmeans[2,:]', "r--", linewidth=3)
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y2, x2nf],["Nonlinear Model Measured","Nonlinear Model No Switch"], loc="best")
xlim([0, tend])
rc("font",size=22)
