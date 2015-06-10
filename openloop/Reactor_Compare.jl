# Compare the transient dynamics of two reactors

include("./params.jl") # load all the parameters and modules

initial_states = [0.5, 400]

xs1 = zeros(2, N)
xs2 = zeros(2, N)

us = ones(N)*0.0
xs1[:,1] = initial_states
xs2[:,1] = initial_states
# Loop through the rest of time
for t=2:N
  if ts[t] < 0.5
    xs1[:, t] = Reactor.run_reactor(xs1[:, t-1], us[t-1], h, cstr_model) # actual plant
    xs2[:, t] = Reactor.run_reactor(xs2[:, t-1], us[t-1], h, cstr_model) # actual plant
  else
    xs1[:, t] = Reactor.run_reactor(xs1[:, t-1], us[t-1], h, cstr_model) # actual plant
    xs2[:, t] = Reactor.run_reactor(xs2[:, t-1], us[t-1], h, cstr_model_broken) # actual plant
  end
end

rc("font", family="serif", size=24)
skip = 50
# figure(1) #
# x1, = plot(xs1[1,:][:], xs1[2,:][:], "k", linewidth=3)
# x2, = plot(xs2[1,:][:], xs2[2,:][:], "r--", linewidth=3)
# ylabel("Temperature [K]")
# xlabel(L"Concentration [kmol.m$^{-3}$]")

figure(2) # Plot filtered results
subplot(2,1,1)
x1, = plot(ts, xs1[1,:]', "k", linewidth=3)
# x2, = plot(ts, xs2[1,:]', "r--", linewidth=3)
ylabel(L"Concentration [kmol.m$^{-3}$]")
xlim([0, tend])
subplot(2,1,2)
x1, = plot(ts, xs1[2,:]', "k", linewidth=3)
# x2, = plot(ts, xs2[2,:]', "r--", linewidth=3)
ylabel("Temperature [K]")
xlabel("Time [min]")
xlim([0, tend])

# legend([x1, x2], [string("CSTR: ", cstr_model.V), string("CSTR: ", cstr_model_broken.V)])
