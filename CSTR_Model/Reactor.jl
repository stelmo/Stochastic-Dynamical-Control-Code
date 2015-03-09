# Run some simulations
# Implement the augmented switching dynamical system
using PyPlot
using Reactor_functions
cd("..\\Linear_Latent_Dynamical_Models")
using LLDS_functions
cd("..\\CSTR_Model")
# Add a definition for convert to make our lives easier!
# But be careful now!
function Base.convert(::Type{Float64}, x::Array{Float64, 1})
  return x[1]
end

# Specify the nonlinear model
cstr_model = begin
  V = 0.1; #m3
  R = 8.314; #kJ/kmol.K
  CA0 = 1.0; #kmol/m3
  TA0 = 310.0; #K
  dH = -4.78e4; #kJ/kmol
  k0 = 72.0e9; #1/min
  E = 8.314e4; #kJ/kmol
  Cp = 0.239; #kJ/kgK
  rho = 1000.0; #kg/m3
  F = 100e-3; #m3/min
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

init_state = [0.57; 395] # initial state
h = 0.001 # time discretisation
tend = 10.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)


initial_states = [0.57, 395]

us = zeros(N)
xs[:,1] = initial_states
# Loop through the rest of time
for t=2:N

  xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], us[t-1], h, cstr_model) # actual plant

  if ts[t] > 0.25
    us[t] = 100.
  end
end

skip = 50
figure(1) # Kalman Filter Demonstration
x1, = plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
ylabel("Temperature [K]")
xlabel(L"Concentration [kmol.m$^{-3}$]")

figure(2) # Plot filtered results
subplot(2,1,1)
x1, = plot(ts, xs[1,:]', "k", linewidth=3)
ylabel(L"Concentration [kmol.m$^{-3}$]")
xlim([0, tend])
subplot(2,1,2)
x2, = plot(ts, xs[2,:]', "k", linewidth=3)
ylabel("Temperature [K]")
xlabel("Time [min]")
xlim([0, tend])
rc("font",size=22)
