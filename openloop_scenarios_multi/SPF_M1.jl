# Inference using two nonlinear models measuring only temperature

tend = 200
include("openloop_params.jl") # load all the parameters and modules

init_state = [0.55, 450]

A = [0.9 0.1;
     0.1 0.9]
# A = [0.5 0.5;0.5 0.5]
fun1(x,u,w) = Reactor.run_reactor(x, u, h, cstr_model) + w
fun2(x,u,w) = Reactor.run_reactor(x, u, h, cstr_model_broken) + w
gs(x) = C1*x
F = [fun1, fun2]
G = [gs, gs]
numSwitches = 2

ydists = [MvNormal(R1);MvNormal(R1)]
xdists = [MvNormal(Q); MvNormal(Q)]
cstr_filter = SPF.Model(F, G, A, xdists, ydists)

nP = 500
xdist = MvNormal(init_state, init_state_covar)
sdist = Categorical([0.9, 0.1])
particles = SPF.init_SPF(xdist, sdist, nP, 2)

switchtrack = zeros(2, N)
maxtrack = zeros(numSwitches, N)
smoothedtrack = zeros(numSwitches, N)

state_noise_dist = MvNormal(Q)
meas_noise_dist = MvNormal(R1)

xs[:,1] = init_state
xsnofix[:,1] = init_state

ys1[1] = C1*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant
SPF.init_filter!(particles, 0.0, ys1[1], cstr_filter)

for k=1:numSwitches
  switchtrack[k, 1] = sum(particles.w[find((x)->x==k, particles.s)])
end
maxtrack[:, 1] = SPF.getMaxTrack(particles, numSwitches)
smoothedtrack[:, 1] = RBPF.smoothedTrack(numSwitches, switchtrack, 1, 10)

spfmeans[:,1], spfcovars[:,:,1] = SPF.getStats(particles)

# Loop through the rest of time
tic()
for t=2:N
  random_element = rand(state_noise_dist)
  if ts[t] < 50.0
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
  else
    xs[:, t] = Reactor.run_reactor(xs[:, t-1], us[t-1], h, cstr_model_broken) + random_element
    xsnofix[:, t] = Reactor.run_reactor(xsnofix[:, t-1], us[t-1], h, cstr_model) + random_element # actual plant
  end

  ys1[t] = C1*xs[:, t] + rand(meas_noise_dist) # measured from actual plant
  SPF.filter!(particles, us[t-1], ys1[t], cstr_filter)
  spfmeans[:,t], spfcovars[:,:,t] = SPF.getStats(particles)

  for k=1:2
    switchtrack[k, t] = sum(particles.w[find((x)->x==k, particles.s)])
  end
  maxtrack[:, t] = SPF.getMaxTrack(particles, numSwitches)
  smoothedtrack[:, t] = RBPF.smoothedTrack(numSwitches, switchtrack, t, 10)

end
toc()
# Plot results
Results.plotSwitchSelection(numSwitches, switchtrack, ts, true)

Results.plotSwitchSelection(numSwitches, maxtrack, ts, false)

Results.plotSwitchSelection(numSwitches, smoothedtrack, ts, false)

Results.plotTrackingBreak(ts, xs, xsnofix, ys1, spfmeans, 1)
