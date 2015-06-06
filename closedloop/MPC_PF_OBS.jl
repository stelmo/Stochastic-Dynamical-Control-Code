# Control using two nonlinear models and measuring both states
# NOTE: remember to adjust the model noise parameter

tend = 150
include("params.jl") # load all the parameters and modules

init_state = [0.5; 450] # initial state

# Set up the PF
f(x, u, w) = Reactor.run_reactor(x, u, h, cstr_model) + w # transistion function
g(x) = C2*x # measurement function
pf_cstr = PF.Model(f,g) # PF object

# Initialise the PF
nP = 100 # number of particles.
prior_dist = MvNormal(init_state, init_state_covar) # prior distribution
particles_pf = PF.init_PF(prior_dist, nP, 2) # initialise the particles

xdist = MvNormal(init_state, init_state_covar)

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
xsnofix[:, 1] = init_state
mnoise = rand(meas_noise_dist)
ys2[:, 1] = C2*xs[:, 1] + mnoise # measured from actual plant
ys2nofix[:, 1] = C2*xsnofix[:, 1] + mnoise # measured from actual plant

PF.init_filter!(particles_pf, 0.0, ys2nofix[:, 1], meas_noise_dist, pf_cstr)
pfmeans[:,1], pfcovars[:,:,1] = PF.getStats(particles_pf)

usnofix[1] = MPC.mpc_mean(pfmeans[:, 1] - lin_models[1].b, horizon, lin_models[1].A, lin_models[1].B, lin_models[1].b, aline, bline, cline, QQ, RR, setpoint[1]-lin_models[1].b[1], usps[1], 15000.0, 1000.0, false)# get the controller input

# Loop through the rest of time
tic()
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

  PF.filter!(particles_pf, usnofix[t-1], ys2nofix[:, t], state_noise_dist, meas_noise_dist, pf_cstr)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles_pf)

  # Controller Input
  if t%10 == 0
    usnofix[t] = MPC.mpc_mean(pfmeans[:, t] - lin_models[1].b, horizon, lin_models[1].A, lin_models[1].B, lin_models[1].b, aline, bline, cline, QQ, RR, setpoint[1]-lin_models[1].b[1], usps[1], 15000.0, 1000.0, false)# get the controller input
  else
    usnofix[t] = usnofix[t-1]
  end
end
toc()
# Plot results
Results.plotTracking(ts, xs, ys2, pfmeans, usnofix, 2, ysp+b[1])
