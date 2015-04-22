# Control using two nonlinear models and measuring both states
# NOTE: remember to adjust the model noise parameter

include("../params.jl") # load all the parameters and modules

initial_states = [0.5; 450] # initial state

# Setup Switching Particle Filter
A = [0.9 0.1;0.1 0.9]
# A = [0.5 0.5;0.5 0.5]
fun1(x,u,w) = Reactor.run_reactor(x, u, h, cstr_model) + w
fun2(x,u,w) = Reactor.run_reactor(x, u, h, cstr_model_broken) + w
gs(x) = C2*x
F = [fun1, fun2]
G = [gs, gs]
numSwitches = 2

ydists = [MvNormal(R2); MvNormal(R2)]
xdists = [MvNormal(Q); MvNormal(Q)]
cstr_filter = SPF.Model(F, G, A, xdists, ydists)

nP = 1000 # number of particles
initial_covar = eye(2)
initial_covar[1] = 1e-3
initial_covar[4] = 4.0
xdist = MvNormal(initial_states, initial_covar)
sdist = Categorical([0.9, 0.1])
particles = SPF.init_SPF(xdist, sdist, nP, 2)

spfmeans = zeros(2, N)
spfcovars = zeros(2,2, N)

switchtrack = zeros(numSwitches, N)
maxtrack = zeros(numSwitches, N)
smoothedtrack = zeros(numSwitches, N)

state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(R2)

# Setup control (use linear control)
linsystems = Reactor.getNominalLinearSystems(h, cstr_model)
linsystems_broken = Reactor.getNominalLinearSystems(h, cstr_model_broken)
opoint = 2 # the specific linear model we will use

lin_models = Array(RBPF.Model, 2)
lin_models[1] = RBPF.Model(linsystems[opoint].A, linsystems[opoint].B, linsystems[opoint].b, C2, Q, R2)
lin_models[2] = RBPF.Model(linsystems_broken[opoint].A, linsystems_broken[opoint].B, linsystems_broken[opoint].b, C2, Q, R2)

H = [1.0 0.0]
controllers = Array(LQR.controller, 2)
for k=1:2
  ysp = 0.4 - lin_models[k].b[1] # set point is set here
  x_off, u_off = LQR.offset(lin_models[k].A, lin_models[k].B, C2, H, ysp)
  K = LQR.lqr(lin_models[k].A, lin_models[k].B, QQ, RR)
  controllers[k] = LQR.controller(K, x_off, u_off)
end

# Setup simulation
xs[:,1] = initial_states
xsnofix[:,1] = initial_states
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant
SPF.init_filter!(particles, 0.0, ys2[:, 1], cstr_filter)

for k=1:2
  switchtrack[k, 1] = sum(particles.w[find((x)->x==k, particles.s)])
end
maxtrack[:, 1] = SPF.getMaxTrack(particles, numSwitches)
smoothedtrack[:, 1] = RBPF.smoothedTrack(numSwitches, switchtrack, 1, 10)

spfmeans[:,1], spfcovars[:,:,1] = SPF.getStats(particles)

# Controller Input
ind = indmax(smoothedtrack[:, 1]) # use this model and controller
us[1] = -controllers[ind].K*(spfmeans[:, 1] - lin_models[ind].b - controllers[ind].x_off) + controllers[ind].u_off # controller action

# Loop through the rest of time
for t=2:N

  random_element = rand(state_noise_dist)
  if ts[t] < 50 # break here
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
  else
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model_broken) + random_element
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
  end

  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  SPF.filter!(particles, us[t-1], ys2[:, t], cstr_filter)
  spfmeans[:,t], spfcovars[:,:,t] = SPF.getStats(particles)

  for k=1:2
    switchtrack[k, t] = sum(particles.w[find((x)->x==k, particles.s)])
  end
  maxtrack[:, t] = SPF.getMaxTrack(particles, numSwitches)
  smoothedtrack[:, t] = RBPF.smoothedTrack(numSwitches, switchtrack, t, 10)

  # Controller Input
  ind = indmax(smoothedtrack[:, t]) # use this model and controller
  us[t] = -controllers[ind].K*(spfmeans[:, t] - lin_models[ind].b - controllers[ind].x_off) + controllers[ind].u_off # controller action

end

# Plot results
rc("font", family="serif", size=24)

figure(1)
axes = Array(Any, numSwitches)
im = 0
width = 500
for k=1:numSwitches
  ax = subplot(numSwitches, 1, k)
  axes[k] = ax
  im = imshow(repeat(smoothedtrack[k,:], outer=[width, 1]), cmap="cubehelix",vmin=0.0, vmax=1.0, interpolation="nearest", aspect="auto")
  tick_params(axis="y", which="both",left="off",right="off", labelleft = "off")
  tick_params(axis="x", which="both",bottom="off", labelbottom = "off")
  ylabel(string("S::",k))
end
tick_params(axis="x", labelbottom = "on")
xticks([1:int(length(ts)/10.0):length(ts)], ts[1:int(length(ts)/10.0):end])
xlabel("Time [min]")


skip = int(length(ts)/20)
skipm = int(length(ts)/20)
figure(3) # Plot filtered results
subplot(3,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
x1nf, = plot(ts, xsnofix[1,:]', "g--", linewidth=3)
y2, = plot(ts[1:skipm:end], ys2[1, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k1, = plot(ts, spfmeans[1,:]', "r--", linewidth=3)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1, k1],["Nonlinear Model","Filtered Mean"], loc="best")
xlim([0, tend])
subplot(3,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
x2nf, = plot(ts, xsnofix[2,:]', "g--", linewidth=3)
y2, = plot(ts[1:skipm:end], ys2[2, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts, spfmeans[2,:]', "r--", linewidth=3)
ylabel("Temperature [K]")
legend([y2, x2nf],["Nonlinear Model Measured","Nonlinear Model No Switch"], loc="best")
xlim([0, tend])
subplot(3,1,3)
plot(ts, us)
ylabel("Controller Input")
xlabel("Time [min]")
xlim([0, tend])
