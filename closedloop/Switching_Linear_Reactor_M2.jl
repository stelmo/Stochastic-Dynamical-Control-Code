# Control using multiple linear models and measuring both concentration and temperature

include("../params.jl") # load all the parameters and modules

init_state = [0.5; 450] # initial state

# Get the three linear models about the nominal operating points
linsystems = Reactor.getNominalLinearSystems(h, cstr_model)

# Setup the RBPF
models, A = RBPF.setup_RBPF(linsystems, C2, Q, R2)
numModels = length(models) # number of linear models (will be 3)
nP = 1000 # number of particles
rbpfmeans = zeros(2, N)
rbpfcovars = zeros(2,2,N)

initial_covar = eye(2) # prior covariance (init_state is the prior mean)
initial_covar[1] = 1e-3
initial_covar[4] = 4.0
sguess =  RBPF.getInitialSwitches(init_state, linsystems) # prior switch distribution
particles = RBPF.init_RBPF(Categorical(sguess), init_state, initial_covar, 2, nP)

maxtrack = zeros(length(linsystems), N) # keep track of the most likely model
switchtrack = zeros(length(linsystems), N) # keep track of the model/switch distribution
smoothedtrack = zeros(length(linsystems), N)

# Setup the controllers
H = [1.0 0.0]
controllers = Array(LQR.controller, length(models))
for k=1:length(models)
  ysp = 0.48 - models[k].b[1] # set point is set here
  x_off, u_off = LQR.offset(models[k].A,models[k].B, C2, H, ysp)
  K = LQR.lqr(models[k].A, models[k].B, QQ, RR)
  controllers[k] = LQR.controller(K, x_off, u_off)
end

state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(R2)
xs[:,1] = init_state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant

RBPF.init_filter!(particles, 0.0, ys2[:, 1], models)
rbpfmeans[:,1], rbpfcovars[:,:, 1] = RBPF.getMLStats(particles)

for k=1:numModels
  switchtrack[k, 1] = sum(particles.ws[find((x)->x==k, particles.ss)])
end
maxtrack[:, 1] = RBPF.getMaxTrack(particles, numModels)
smoothedtrack[:, 1] = RBPF.smoothedTrack(numModels, switchtrack, 1, 10)

# Controller Input
# insert (hopefully)!

#Loop through the rest of time
for t=2:N

  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  RBPF.filter!(particles, us[t-1], ys2[:, t], models, A)
  rbpfmeans[:,t], rbpfcovars[:,:, t] = RBPF.getMLStats(particles)

  for k=1:numModels
    switchtrack[k, t] = sum(particles.ws[find((x)->x==k, particles.ss)])
  end
  maxtrack[:, t] = RBPF.getMaxTrack(particles, numModels)
  smoothedtrack[:, t] = RBPF.smoothedTrack(numModels, switchtrack, t, 5)

  # Controller Input
  # ind = indmax(maxtrack[:, t-1]) # use this model and controller
  # us[t-1] = -controllers[ind].K*(rbpfmeans[:, t-1] - models[ind].b - controllers[ind].x_off) + controllers[ind].u_off # controller action
  # if us[t-1] > 6000
  #   us[t-1] = 6000
  # end
  # if us[t-1] < -6000
  #   us[t-1] = -6000
  # end

end


# Plot results
rc("font", family="serif", size=24)

figure(1) # Model and state space
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

figure(2) # Model selection
axes = Array(Any, length(linsystems))
im = 0
width = 500
for k=1:length(linsystems)
  ax = subplot(length(linsystems), 1, k)
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

figure(3) # Model selection
axes = Array(Any, length(linsystems))
im = 0
width = 500
for k=1:length(linsystems)
  ax = subplot(length(linsystems), 1, k)
  axes[k] = ax
  im = imshow(repeat(smoothedtrack[k,:], outer=[width, 1]), cmap="cubehelix",vmin=0.0, vmax=1.0, interpolation="nearest", aspect="auto")
  tick_params(axis="y", which="both",left="off",right="off", labelleft = "off")
  tick_params(axis="x", which="both",bottom="off", labelbottom = "off")
  ylabel(string("S::",k))
end
tick_params(axis="x", labelbottom = "on")
xticks([1:int(length(ts)/10.0):length(ts)], ts[1:int(length(ts)/10.0):end])
# colorbar(im, ax=axes)
xlabel("Time [min]")

figure(4) # Tracking
skip = int(length(ts)/75)
skipm = skip
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
y2, = plot(ts[1:skipm:end], ys2[1, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k1, = plot(ts[1:skip:end], rbpfmeans[1, 1:skip:end]', "r--",linewidth=3)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1],["Nonlinear Model"], loc="best")
xlim([0, tend])
ylim([0, 1])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
y2, = plot(ts[1:skipm:end], ys2[2, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts[1:skip:end], rbpfmeans[2, 1:skip:end]', "r--",linewidth=3)
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y2, k2],["Nonlinear Model Measured", "Filtered Mean Estimate"], loc="best")
xlim([0, tend])
