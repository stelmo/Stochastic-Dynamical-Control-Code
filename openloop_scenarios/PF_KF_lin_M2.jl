# Inference using a linear model in the PF compared to the KF (which implicitly uses a linear model)
# NOTE: don't run for too long because the linear model is unstable!

tend = 20
include("openloop_params.jl") # load all the parameters and modules

init_state = [0.50, 400]

# Specify the linear model
linsystems = Reactor.getNominalLinearSystems(h, cstr_model)
opoint = 2
A = linsystems[opoint].A
B = linsystems[opoint].B
b = linsystems[opoint].b

lin_cstr = LLDS.llds(A, B, C2, Q, R2) # KF object

# Setup the PF
f(x, u, w) = A*x + B*u + w
g(x) = C2*x # state observation
cstr_pf = PF.Model(f,g)
nP = 500 #number of particles.

# Initialise the PFs
init_pf_dist = MvNormal(init_state-b, init_state_covar) # prior distribution
particles = PF.init_PF(init_pf_dist, nP, 2) # initialise the particles

state_noise_dist = MvNormal(Q) # state distribution
meas_noise_dist = MvNormal(R2) # measurement distribution

pfmeans = zeros(2, N)
pfcovars = zeros(2,2, N)
kfmeans = zeros(2, N)
kfcovars = zeros(2, 2, N)

# Time step 1
xs[:,1] = init_state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys2[:, 1]-b, meas_noise_dist, cstr_pf)
pfmeans[:,1], pfcovars[:,:,1] = PF.getStats(particles)
kfmeans[:, 1], kfcovars[:, :, 1] = LLDS.init_filter(init_state-b, init_state_covar, ys2[:, 1]-b, lin_cstr)

# Loop through the rest of time
for t=2:N
  xs[:, t] = A*(xs[:, t-1]-b) + B*us[t-1] + b + rand(state_noise_dist)# actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  PF.filter!(particles, us[t-1], ys2[:, t]-b, state_noise_dist, meas_noise_dist, cstr_pf)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles)
  kfmeans[:, t], kfcovars[:,:, t] = LLDS.step_filter(kfmeans[:, t-1], kfcovars[:,:, t-1], us[t-1], ys2[:, t]-b, lin_cstr)
end

pfmeans = pfmeans .+ b
kfmeans = kfmeans .+ b

# Plot Results
Results.plotEllipseComp(pfmeans, pfcovars, kfmeans, kfcovars, xs, ts)

Results.plotTrackingTwoFilters(ts, xs, ys2, pfmeans, kfmeans)

println("For the Kalman Filter:")
avediff = Results.calcError(xs, kfmeans)
println("For the Particle Filter:")
avediff = Results.calcError(xs, pfmeans)
