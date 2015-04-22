# Control using two nonlinear models and measuring both states

include("../params.jl") # load all the parameters and modules

initial_states = [0.5; 400] # initial state


# Setup Switching Particle Filter
A = [0.9 0.1;0.1 0.9]
# A = [0.5 0.5;0.5 0.5]
fun1(x,u,w) = Reactor.run_reactor(x, u, h, cstr_model)
fun2(x,u,w) = Reactor.run_reactor(x, u, h, cstr_model_broken)
gs(x) = newC*x
F = [fun1, fun2]
G = [gs, gs]
numSwitches = 2

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
maxtrack = zeros(numSwitches, N)

state_dist = MvNormal(Q)

measurements = MvNormal(R)
xs[:,1] = initial_states
xsnofix[:,1] = initial_states

ys[:, 1] = newC*xs[:, 1] + rand(measurements) # measured from actual plant
SPF.init_filter!(particles, 0.0, ys[:, 1], cstr_filter)

for k=1:2
  switchtrack[k, 1] = sum(particles.w[find((x)->x==k, particles.s)])
end
maxtrack[:, 1] = SPF.getMaxTrack(particles, numSwitches)

fmeans[:,1], fcovars[:,:,1] = SPF.getStats(particles)

# Loop through the rest of time
for t=2:N

  random_element = rand(state_dist)
  if ts[t] < 5.0
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
  else
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model_broken) + random_element
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
  end

  ys[:, t] = newC*xs[:, t] + rand(measurements) # measured from actual plant
  SPF.filter!(particles, us[t-1], ys[:, t], cstr_filter)
  fmeans[:,t], fcovars[:,:,t] = SPF.getStats(particles)

  for k=1:2
    switchtrack[k, t] = sum(particles.w[find((x)->x==k, particles.s)])
  end
  maxtrack[:, t] = SPF.getMaxTrack(particles, numSwitches)

end

rc("font", family="serif", size=24)

figure(1)
axes = Array(Any, numSwitches)
im = 0
width = 500
for k=1:numSwitches
  ax = subplot(numSwitches, 1, k)
  axes[k] = ax
  im = imshow(repeat(maxtrack[k,:], outer=[width, 1]), cmap="cubehelix",vmin=0.0, vmax=1.0, interpolation="nearest", aspect="auto")
  tick_params(axis="y", which="both",left="off",right="off", labelleft = "off")
  tick_params(axis="x", which="both",bottom="off", labelbottom = "off")
  ylabel(string("S::",k))
end
tick_params(axis="x", labelbottom = "on")
xticks([1:int(length(ts)/10.0):length(ts)], ts[1:int(length(ts)/10.0):end])
# colorbar(im, ax=axes)
xlabel("Time [min]")

# figure(1)
# maxswitch = maximum(switchtrack)
# axes = Array(Any, 2)
# im = 0
# width = 500
# for k=1:2
#   ax = subplot(2, 1, k)
#   axes[k] = ax
#   im = imshow(repeat(switchtrack[k,:], outer=[width, 1]), cmap="cubehelix",vmin=0.0, vmax=maxswitch, interpolation="nearest", aspect="auto")
#   tick_params(axis="y", which="both",left="off",right="off", labelleft = "off")
#   tick_params(axis="x", which="both",bottom="off", labelbottom = "off")
#   ylabel(string("S::",k))
# end
# tick_params(axis="x", labelbottom = "on")
# xticks([1:int(length(ts)/10.0):length(ts)], ts[1:int(length(ts)/10.0):end])
# colorbar(im, ax=axes)
# xlabel("Time [min]")


skip = 50
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
rc("font",size=22)
