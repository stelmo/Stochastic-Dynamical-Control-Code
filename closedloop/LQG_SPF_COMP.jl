# Control using two nonlinear models and measuring both states
# NOTE: remember to adjust the model noise parameter

include("../params.jl") # load all the parameters and modules

init_state = [0.5; 450] # initial state

# Setup Switching Particle Filter
A = [0.9 0.1;0.1 0.9]
# A = [0.5 0.5;0.5 0.5]
fun1(x,u,w) = Reactor.run_reactor(x, u, h, cstr_model) + w
fun2(x,u,w) = Reactor.run_reactor(x, u, h, cstr_model_broken) + w
gs(x) = C2*x
F = [fun1, fun2]
G = [gs, gs]
numSwitches = 2

# Set up the PF
f(x, u, w) = Reactor.run_reactor(x, u, h, cstr_model) + w # transistion function
g(x) = C2*x # measurement function
pf_cstr = PF.Model(f,g) # PF object
# Initialise the PF
nP = 100 # number of particles.
prior_dist = MvNormal(init_state, init_state_covar) # prior distribution
particles_pf = PF.init_PF(prior_dist, nP, 2) # initialise the particles


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

H = [1.0 0.0]
controllers = Array(LQR.controller, 2)
for k=1:2
  ysp = 0.48 - lin_models[k].b[1] # set point is set here
  x_off, u_off = LQR.offset(lin_models[k].A, lin_models[k].B, C2, H, ysp)
  K = LQR.lqr(lin_models[k].A, lin_models[k].B, QQ, RR)
  controllers[k] = LQR.controller(K, x_off, u_off)
end

# Setup simulation
xs[:, 1] = init_state
xsnofix[:, 1] = init_state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant
ys2nofix[:, 1] = C2*xsnofix[:, 1] + rand(meas_noise_dist) # measured from actual plant

SPF.init_filter!(particles, 0.0, ys2[:, 1], cstr_filter)

PF.init_filter!(particles_pf, 0.0, ys2nofix[:, 1], meas_noise_dist, pf_cstr)
pfmeans[:,1], pfcovars[:,:,1] = PF.getStats(particles_pf)

for k=1:2
  switchtrack[k, 1] = sum(particles.w[find((x)->x==k, particles.s)])
end
maxtrack[:, 1] = SPF.getMaxTrack(particles, numSwitches)
smoothedtrack[:, 1] = RBPF.smoothedTrack(numSwitches, switchtrack, 1, 10)

spfmeans[:,1], spfcovars[:,:,1] = SPF.getStats(particles)

# Controller Input
ind = indmax(smoothedtrack[:, 1]) # use this model and controller
us[1] = -controllers[ind].K*(spfmeans[:, 1] - lin_models[ind].b - controllers[ind].x_off) + controllers[ind].u_off # controller action
usnofix[1] = -controllers[1].K*(pfmeans[:, 1] - lin_models[1].b - controllers[1].x_off) + controllers[1].u_off

# Loop through the rest of time
for t=2:N

  random_element = rand(state_noise_dist)
  if ts[t] < 50 # break here
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], usnofix[t-1], h, cstr_model) + random_element # actual plant
  else
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model_broken) + random_element
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], usnofix[t-1], h, cstr_model_broken) + random_element # actual plant
  end

  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  ys2nofix[:, t] = C2*xsnofix[:, t] + rand(meas_noise_dist) # measured from actual plant

  SPF.filter!(particles, us[t-1], ys2[:, t], cstr_filter)
  spfmeans[:,t], spfcovars[:,:,t] = SPF.getStats(particles)

  PF.filter!(particles_pf, usnofix[t-1], ys2nofix[:, t], state_noise_dist, meas_noise_dist, pf_cstr)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles_pf)


  for k=1:2
    switchtrack[k, t] = sum(particles.w[find((x)->x==k, particles.s)])
  end
  maxtrack[:, t] = SPF.getMaxTrack(particles, numSwitches)
  smoothedtrack[:, t] = RBPF.smoothedTrack(numSwitches, switchtrack, t, 20)

  # Controller Input
  ind = indmax(smoothedtrack[:, t]) # use this model and controller
  us[t] = -controllers[ind].K*(spfmeans[:, t] - lin_models[ind].b - controllers[ind].x_off) + controllers[ind].u_off # controller action
  usnofix[t] = -controllers[1].K*(pfmeans[:, t] - lin_models[1].b - controllers[1].x_off) + controllers[1].u_off # controller action

end

# Plot results
Results.plotSwitchSelection(numSwitches, smoothedtrack, ts, false)

Results.plotTrackingComparison(ts, xs, us, xsnofix, usnofix)
