# Inference using two nonlinear models measuring only temperature

include("../params.jl") # load all the parameters and modules

init_state = [0.50, 400]

A = [0.9 0.1;0.1 0.9]
# A = [0.5 0.5;0.5 0.5]
fun1(x,u,w) = Reactor.run_reactor(x, u, h, cstr_model) + w
fun2(x,u,w) = Reactor.run_reactor(x, u, h, cstr_model_broken) + w
gs(x) = C2*x
F = [fun1, fun2]
G = [gs, gs]

ydists = [MvNormal(R2);MvNormal(R2)]
xdists = [MvNormal(Q); MvNormal(Q)]
cstr_filter = SPF.Model(F, G, A, xdists, ydists)

nP = 1500
xdist = MvNormal(init_state, init_state_covar)
sdist = Categorical([0.5, 0.5])
particles = SPF.init_SPF(xdist, sdist, nP, 2)

switchtrack = zeros(2, N)
maxtrack = zeros(2, N)
smoothedtrack = zeros(2, N)

state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(R2)

xs[:,1] = init_state
xsnofix[:,1] = init_state

ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant
SPF.init_filter!(particles, 0.0, ys2[:, 1], cstr_filter)

for k=1:2
  switchtrack[k, 1] = sum(particles.w[find((x)->x==k, particles.s)])
end
maxtrack[:, 1] = SPF.getMaxTrack(particles, numSwitches)
smoothedtrack[:, 1] = RBPF.smoothedTrack(numSwitches, switchtrack, 1, 10)

spfmeans[:,1], spfcovars[:,:,1] = SPF.getStats(particles)
# Loop through the rest of time
for t=2:N
  random_element = rand(state_noise_dist)
  if ts[t] < 50.0
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


end

# Plot results
Results.plotSwitchSelection(numSwitches, switchtrack, ts, true)

Results.plotSwitchSelection(numSwitches, maxtrack, ts, false)

Results.plotSwitchSelection(numSwitches, smoothedtrack, ts, false)

Results.plotTrackingBreak(ts, xs, xsnofix, ys2, spfmeans, 2)
