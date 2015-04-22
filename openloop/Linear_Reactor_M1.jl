# Inference using one linear model measuring only temperature

include("../params.jl") # load all the parameters and modules

init_state = [0.50, 400]

# Specify the linear model
linsystems = Reactor.getNominalLinearSystems(h, cstr_model)
opoint = 2 # which nominal model to use
A = linsystems[opoint].A
B = linsystems[opoint].B
b = linsystems[opoint].b

lin_cstr = LLDS.llds(A, B, C1, Q, R1) # KF object

# Plant initialisation
xs[:,1] = init_state
linxs[:, 1] = init_state - b

# Simulate plant
state_noise_dist = MvNormal(Q)
meas_noise_dist = Normal(0.0, sqrt(R1))
ys1[1] = C1*xs[:, 1] + rand(meas_noise_dist) # measure from actual plant

# Filter setup
kfmeans = zeros(2, N)
kfcovars = zeros(2,2, N)
init_mean = init_state - b
init_covar = eye(2) # vague
init_covar[1] = 1e-3
init_covar[4] = 10.

# First time step
kfmeans[:, 1], kfcovars[:,:, 1] = LLDS.init_filter(init_mean, init_covar, ys1[1]-b[2], lin_cstr)

for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist)# actual plant
  ys1[t] = C1*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  linxs[:, t], temp = LLDS.step(linxs[:, t-1], us[t-1], lin_cstr)
  kfmeans[:, t], kfcovars[:,:, t] = LLDS.step_filter(kfmeans[:, t-1], kfcovars[:,:, t-1], us[t-1], ys1[t] - b[2], lin_cstr)
end

linxs = linxs .+ b
kfmeans = kfmeans .+ b

# Plot results

rc("font", family="serif", size=24)

skip = int(length(ts)/20)
figure(1)
x1, = plot(xs[1,:][:], xs[2,:][:], "k",linewidth=3)
f1, = plot(kfmeans[1, 1:skip:end][:], kfmeans[2, 1:skip:end][:], "bx", markersize=5, markeredgewidth = 2)
b1 = 0.0
for k=1:skip:N
  p1, p2 = Ellipse.ellipse(kfmeans[:,k], kfcovars[:,:, k])
  b1, = plot(p1, p2, "b")
end
plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
ylabel("Temperature [K]")
xlabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1,f1, b1],["Nonlinear Model","Kalman Filter Mean", L"Kalman Filter $1\sigma$-Ellipse"], loc="best")


skipm = int(length(ts)/20)
figure(2) # Filtering
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
linx1, = plot(ts, linxs[1,:]', "r--", linewidth=3)
k1, = plot(ts[1:skip:end], kfmeans[1, 1:skip:end]', "bx", markersize=5, markeredgewidth = 2)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([x1,linx1],["Nonlinear Model","Linear Model"], loc="best")
xlim([0, tend])
ylim([0, 1.5])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
linx2, = plot(ts, linxs[2,:]', "r--", linewidth=3)
y1, = plot(ts[1:skipm:end], ys1[1:skipm:end], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts[1:skip:end], kfmeans[2, 1:skip:end]', "bx", markersize=5, markeredgewidth = 2)
ylabel("Temperature [K]")
xlabel("Time [min]")
legend([y1, k2],["Nonlinear Model Measured", "Filtered Mean Estimate"], loc="best")
xlim([0, tend])
ylim([minimum(xs[2,:]), maximum(xs[2,:])])
