# Inference using a linear model in the PF compared to the KF (which implicitly uses a linear model)
# NOTE: don't run for too long because the linear model is unstable!

include("../params.jl") # load all the parameters and modules

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
init_state_covar = eye(2)*1e-3 # initial covariance
init_state_covar[4] = 4.0
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
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist)# actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  PF.filter!(particles, us[t-1], ys2[:, t]-b, state_noise_dist, meas_noise_dist, cstr_pf)
  pfmeans[:,t], pfcovars[:,:,t] = PF.getStats(particles)
  kfmeans[:, t], kfcovars[:,:, t] = LLDS.step_filter(kfmeans[:, t-1], kfcovars[:,:, t-1], us[t-1], ys2[:, t]-b, lin_cstr)
end

pfmeans = pfmeans .+ b
kfmeans = kfmeans .+ b

rc("font", family="serif", size=24)

skip = int(length(ts)/20)
figure(1) #  Filter Demonstration
x1, = plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
f1, = plot(pfmeans[1, 1:skip:end][:], pfmeans[2, 1:skip:end][:], "rx", markersize=5, markeredgewidth = 2)
f2, = plot(kfmeans[1, 1:skip:end][:], kfmeans[2, 1:skip:end][:], "bx", markersize=5, markeredgewidth = 2)
b1 = 0.0
b2 = 0.0
for k=1:skip:N
  p1, p2 = Ellipse.ellipse(pfmeans[:,k], pfcovars[:,:, k])
  b1, = plot(p1, p2, "r")

  p3, p4 = Ellipse.ellipse(kfmeans[:,k], kfcovars[:,:, k])
  b2, = plot(p3, p4, "b")
end
plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
ylabel("Temperature [K]")
xlabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1,f1,f2, b1, b2],["Nonlinear Model","Particle Filter Mean","Kalman Filter Mean", L"Particle Filter $1\sigma$-Ellipse", L"Kalman Filter $1\sigma$-Ellipse"], loc="best")

skipm = int(length(ts)/20)
figure(2) # Plot filtered results
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
k1, = plot(ts[1:skip:end], pfmeans[1,1:skip:end]', "rx", markersize=5, markeredgewidth=2)
y2, = plot(ts[1:skipm:end], ys2[1, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k12, = plot(ts[1:skip:end], kfmeans[1, 1:skip:end]', "bx", markersize=5, markeredgewidth=2)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1, k1],["Nonlinear Model","Particle Filter"], loc="best")
xlim([0, tend])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
y2, = plot(ts[1:skipm:end], ys2[2, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts[1:skip:end], pfmeans[2,1:skip:end]', "rx", markersize=5, markeredgewidth=2)
k22, = plot(ts[1:skip:end], kfmeans[2, 1:skip:end]', "bx", markersize=5, markeredgewidth=2)
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y2, k22],["Nonlinear Model Measured", "Kalman Filter"], loc="best")
xlim([0, tend])
