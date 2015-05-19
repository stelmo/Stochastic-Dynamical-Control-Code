# Control using two nonlinear models and measuring both states
# NOTE: remember to adjust the model noise parameter
tend = 50
include("params.jl") # load all the parameters and modules

init_state = [0.5; 450] # initial state

# Setup Switching Particle Filter
A = [0.999 0.001;0.001 0.999]

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
xdist = MvNormal(init_state, init_state_covar)
sdist = Categorical([0.9, 0.1])
particles = SPF.init_SPF(xdist, sdist, nP, 2)

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

setpoint = linsystems[opoint].op
H = [1.0 0.0] # only attempt to control the concentration
usps = zeros(length(lin_models))
for k=1:length(lin_models)
  sp = setpoint[1] - lin_models[k].b[1]
  x_off, usp = LQR.offset(lin_models[k].A, lin_models[k].B, C2, H, sp) # control offset
  usps[k] = usp[1]
end

horizon = 150
# add state constraints
aline = 10. # slope of constraint line ax + by + c = 0
cline = -403.0 # negative of the y axis intercept
bline = 1.0

# Setup simulation
xs[:, 1] = init_state
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
yspfix = setpoint[1] - lin_models[ind].b[1]
us[1] = MPC.mpc_mean(spfmeans[:, 1] - lin_models[ind].b, horizon, lin_models[ind].A, lin_models[ind].B, lin_models[ind].b, aline, bline, cline, QQ, RR, yspfix, usps[ind], 15000.0, 1000.0, false)# get the controller input

# Loop through the rest of time
tic()
for t=2:N

  random_element = rand(state_noise_dist)
  if ts[t] < 50 # break here
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
  else
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model_broken) + random_element
  end

  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant

  SPF.filter!(particles, us[t-1], ys2[:, t], cstr_filter)
  spfmeans[:,t], spfcovars[:,:,t] = SPF.getStats(particles)

  for k=1:2
    switchtrack[k, t] = sum(particles.w[find((x)->x==k, particles.s)])
  end
  maxtrack[:, t] = SPF.getMaxTrack(particles, numSwitches)
  smoothedtrack[:, t] = RBPF.smoothedTrack(numSwitches, switchtrack, t, 40)

  # Controller Input
  ind = indmax(smoothedtrack[:, t]) # use this model and controller
  yspfix = setpoint[1] - lin_models[ind].b[1]
  us[t] = MPC.mpc_mean(spfmeans[:, t] - lin_models[ind].b, horizon, lin_models[ind].A, lin_models[ind].B, lin_models[ind].b, aline, bline, cline, QQ, RR, yspfix, usps[ind], 15000.0,1000.0, false)# get the controller input

end
toc()

# Plot results
Results.plotSwitchSelection(numSwitches, maxtrack, ts, false)

Results.plotSwitchSelection(numSwitches, smoothedtrack, ts, false)

Results.plotTracking(ts, xs, ys2, spfmeans, us, 2, setpoint[1])

Results.plotEllipses(ts, xs, spfmeans, spfcovars, "MPC", [aline, cline], setpoint, true, 4.6052)
