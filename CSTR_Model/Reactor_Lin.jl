# Linearisation Procedure
using PyPlot
import Reactor_functions
reload("Reactor_functions.jl")

# Introduce the reactor
cstr = begin
  V = 5.0 #m3
  R = 8.314 #kJ/kmol.K
  CA0 = 1.0 #kmol/m3
  TA0 = 310.0 #K
  dH = -4.78e4 #kJ/kmol
  k0 = 72.0e7 #1/min
  E = 8.314e4 #kJ/kmol
  Cp = 0.239 #kJ/kgK
  rho = 1000.0 #kg/m3
  F = 100e-3 #m3/min
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

h = 0.001 # time discretisation
tend = 1.0 # end simulation time
ts = [0.0:h:tend]

# Divide state space into sectors: n by m
nX = 3 # rows
nY = 3 # cols
xspace = [0.0, 1.0]
yspace = [250, 550]

linsystems = Reactor_functions.getLinearSystems(nX, nY, xspace, yspace, h, cstr)

# figure(1)
for k=1:(nX*nY)
  initial_states = linsystems[k].op + randn(2).*[0.01, 2.0]

  N = length(ts)
  xs = zeros(2, N)
  linxs = zeros(2, N)
  xs[:,1] = initial_states
  linxs[:,1] = initial_states
  # Loop through the rest of time
  for t=2:N
      xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], 0.0, h, cstr) # actual plant
      linxs[:, t] = linsystems[k].A*linxs[:, t-1] + linsystems[k].B*0.0 + linsystems[k].b
  end

  subplot(nX, nY, k)
  # figure(k)
  x1, = plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
  x2, = plot(linxs[1,:][:], linxs[2,:][:], "r--", linewidth=3)
  plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
  plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
  # xlim([minimum(xs[1,:])-0.1, maximum(xs[1,:])+0.1])
  # xlim([0.0, 1.0])
  # ylim([minimum(xs[2,:])-10.0, maximum(xs[2,:])+10.0])
  # (k == 1 || k== 2) && ylabel("Temperature [K]")
  # xlabel(L"Concentration [kmol.m$^{-3}$]")
end
