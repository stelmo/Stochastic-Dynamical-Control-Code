# Controller using the linear reactor model measuring both concentration and temperature.

include("../params.jl") # load all the parameters and modules

# Get the linear model
linsystems = Reactor.getNominalLinearSystems(h, cstr_model) # cstr_model comes from params.jl
opoint = 2 # the specific operating point we are going to use for control

init_state = linsystems[opoint].op + rand(-1:2:1)*rand(2).*[0.2, 10] # random initial point near operating point

# Set the state space model
A = linsystems[opoint].A
B = linsystems[opoint].B
b = linsystems[opoint].b # offset from the origin

# Set point
# ysp = linsystems[1].op[1] - b[1] # Low concentration
ysp = linsystems[2].op[1] - b[1] # Medium concentration
# ysp = linsystems[3].op[1] - b[1] # High concentration
ysp = 0.7 - b[1] # Medium concentration

# Create the controller
H = [1.0 0.0] # only attempt to control the concentration
x_off, u_off = LQR.offset(A,B,C2,H, ysp) # control offset
K = LQR.lqr(A, B, QQ, RR) # controller

# Set up the KF
kf_cstr = LLDS.llds(A, B, C2, Q, R2) # set up the KF object (measuring both states)
state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(R2)
kfmeans = zeros(2, N)
kfcovars = zeros(2,2, N)

# Filter initialisation
init_mean = init_state - b # state space offset
init_covar = eye(2)
init_covar[1] = 1e-3
init_covar[4] = 4.

# First time step of the simulation
xs[:,1] = init_state # set simulation starting point to the random initial state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measure from actual plant
kfmeans[:, 1], kfcovars[:,:, 1] = LLDS.init_filter(init_mean, init_covar, ys2[:, 1]-b, kf_cstr) # filter
us[1] = -K*(kfmeans[:, 1] - x_off) + u_off # controller action

for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + rand(state_noise_dist) # actual plant
  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist) # measure from actual plant
  kfmeans[:, t], kfcovars[:,:, t] = LLDS.step_filter(kfmeans[:, t-1], kfcovars[:,:, t-1], us[t-1], ys2[:, t]-b, kf_cstr)

  # Compute controller action
  us[t] = -K*(kfmeans[:, t] - x_off) + u_off # controller action
  # if t%10 == 0
  #   us[t-1] = -K*(kfmeans[:, t-1] - x_off) + u_off # controller action
  # else
  #   us[t-1] = us[t-2]
  # end

end
kfmeans = kfmeans .+ b

# Plot the results
rc("font", family="serif", size=24)

skipmeas = int(length(ts)/20)
skipmean = int(length(ts)/20)
figure(1) # Filtering
subplot(3,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
y2, = plot(ts[1:skipmeas:end], ys2[1, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
k1, = plot(ts[1:skipmean:end], kfmeans[1, 1:skipmean:end]', "bx", markersize=5, markeredgewidth = 2)
ylabel(L"Concentration [kmol.m$^{-3}$]")
legend([y2, x1],["Nonlinear Model Measured", "Nonlinear Model"], loc="best")
xlim([0, tend])
ylim([0, 1])

subplot(3,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
y2, = plot(ts[1:skipmeas:end], ys2[2, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
k2, = plot(ts[1:skipmean:end], kfmeans[2, 1:skipmean:end]', "bx", markersize=5, markeredgewidth = 2)
ylabel("Temperature [K]")
legend([k2],["Filtered Mean Estimate"], loc="best")
xlim([0, tend])
ylim([minimum(xs[2,:]), maximum(xs[2,:])])
subplot(3,1,3)
plot(ts, us)
xlim([0, tend])
xlabel("Time [min]")
ylabel("Controller Input")
