# Inference using two nonlinear models measuring temperature vs measuring temperature and concentration

include("./params.jl") # load all the parameters and modules

init_state = [0.50, 400]

A = [0.9 0.1;0.1 0.9]
# A = [0.5 0.5;0.5 0.5]

fun1(x,u,w) = Reactor.run_reactor(x, u, h, cstr_model) + w
fun2(x,u,w) = Reactor.run_reactor(x, u, h, cstr_model_broken) + w
F = [fun1, fun2]

gs1(x) = C1*x
gs2(x) = C2*x
G1 = [gs1, gs1]
G2 = [gs2, gs2]

xdists = [MvNormal(Q), MvNormal(Q)]

ydists1 = [MvNormal(R1), MvNormal(R1)]
cstr_filter1 = SPF.Model(F, G1, A, xdists, ydists1)

ydists2 = [MvNormal(R2), MvNormal(R2)]
cstr_filter2 = SPF.Model(F, G2, A, xdists, ydists2)


nP = 1500
xdist = MvNormal(init_state, init_state_covar)
sdist = Categorical([0.5, 0.5])

particles1 = SPF.init_SPF(xdist, sdist, nP, 2)
spfmeans1 = zeros(2, N)
spfcovars1 = zeros(2, 2, N)

particles2 = SPF.init_SPF(xdist, sdist, nP, 2)
spfmeans2 = zeros(2, N)
spfcovars2 = zeros(2, 2, N)

meas_noise_dist1 = MvNormal(R1)
meas_noise_dist2 = MvNormal(R2)
xs[:, 1] = init_state
xsnofix[:, 1] = init_state

ys1[1] = C1*xs[:, 1] + rand(meas_noise_dist1) # measured from actual plant
SPF.init_filter!(particles1, 0.0, ys1[1], cstr_filter1)
spfmeans1[:,1], spfcovars1[:,:,1] = SPF.getStats(particles1)

ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist2) # measured from actual plant
SPF.init_filter!(particles2, 0.0, ys2[:, 1], cstr_filter2)
spfmeans2[:,1], spfcovars2[:,:,1] = SPF.getStats(particles2)

# Loop through the rest of time
for t=2:N
  temp = rand(xdists[1])
  if ts[t] < 50.0
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + temp # actual plant
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr_model) + temp # actual plant
  else
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model_broken) + temp
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr_model) +temp # actual plant
  end

  ys1[t] = C1*xs[:, t] + rand(meas_noise_dist1) # measured from actual plant
  SPF.filter!(particles1, us[t-1], ys1[t], cstr_filter1)
  spfmeans1[:,t], spfcovars1[:,:,t] = SPF.getStats(particles1)

  ys2[:, t] = C2*xs[:, t] + rand(meas_noise_dist2) # measured from actual plant
  SPF.filter!(particles2, us[t-1], ys2[:, t], cstr_filter2)
  spfmeans2[:,t], spfcovars2[:,:,t] = SPF.getStats(particles2)

end


# Plot Results
Results.plotEllipseComp(spfmeans1, spfcovars1, "M1", spfmeans2, spfcovars2, "M2", xs, ts)
