# Linearisation Procedure
using PyPlot
import Reactor_functions
reload("Reactor_functions.jl")

# Introduce the reactor
cstr = begin
  V = 0.1 #m3
  R = 8.314 #kJ/kmol.K
  CA0 = 1.0 #kmol/m3
  TA0 = 310.0 #K
  dH = -4.78e4 #kJ/kmol
  k0 = 72.0e9 #1/min
  E = 8.314e4 #kJ/kmol
  Cp = 0.239 #kJ/kgK
  rho = 1000.0 #kg/m3
  F = 100e-3 #m3/min
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

h = 0.001 # time discretisation
tend = 5. # end simulation time
ts = [0.0:h:tend]

# Divide state space into sectors: n by m
nX = 2 # rows
nY = 3 # cols
xspace = [0.0, 1.0]
yspace = [250, 650]

linsystems = Reactor_functions.getLinearSystems(nX, nY, xspace, yspace, h, cstr)

# linpoint = [0.4, 450]
# A, B, b = Reactor_functions.linearise(linpoint, h, cstr)
# N = length(ts)
# xs1 = zeros(2, N)
# linxs = zeros(2, N)
#
# initial_states = linpoint + [0.05, 5.0]
#
# us = ones(N)*0.0
# xs1[:,1] = initial_states
# linxs[:,1] = initial_states
# # Loop through the rest of time
# for t=2:N
#     xs1[:, t] = Reactor_functions.run_reactor(xs1[:, t-1], us[t-1], h, cstr) # actual plant
#     linxs[:, t] = A*linxs[:, t-1] + B*us[t-1] + b
# end
#
# skip = 50
# figure(1) # Kalman Filter Demonstration
# x1, = plot(xs1[1,:][:], xs1[2,:][:], "k", linewidth=3)
# x2, = plot(linxs[1,:][:], linxs[2,:][:], "r--", linewidth=3)
# plot(xs1[1,1], xs1[2,1], "ko", markersize=10, markeredgewidth = 4)
# plot(xs1[1,end], xs1[2,end], "kx", markersize=10, markeredgewidth = 4)
# xlim([minimum(xs1[1,:])-0.1, maximum(xs1[1,:])+0.1])
# ylim([minimum(xs1[2,:])-10.0, maximum(xs1[2,:])+10.0])
# ylabel("Temperature [K]")
# xlabel(L"Concentration [kmol.m$^{-3}$]")
