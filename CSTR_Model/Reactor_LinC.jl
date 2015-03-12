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
tend = 50.0 # end simulation time
ts = [0.0:h:tend]

# Divide state space into sectors: n by m
npoints = 1
total_ops = npoints + 3 # nX*nY +3

xspace = [-0.5, 1.0]
yspace = [250, 650]

# linsystems = Reactor_functions.getLinearSystems(nX, nY, xspace, yspace, h, cstr)
linsystems = Reactor_functions.getLinearSystems_randomly(npoints, xspace, yspace, h, cstr)

figure(1)
k=4 # set which operating point to use...
nDD = 2
vars = zeros(2,3)
vars[:,1] = [0.05, 5.0]
vars[:,2] = [0.3, 20.0]
# vars[:,3] = [0.08, 16.0]
x1 = 0
x2 = 0
x3 = 0
for dd=1:nDD # only loop through
  initial_states = linsystems[2].op + abs(randn(2)).*vars[:, dd]

  N = length(ts)
  xs = zeros(2, N)
  linxs = zeros(2, N)
  xs[:,1] = initial_states
  linxs[:,1] = initial_states
  # Loop through the rest of time
  for t=2:N
      xs[:, t] = Reactor_functions.run_reactor(xs[:, t-1], 0.0, h, cstr) # actual plant
      temp = linsystems[k].A*linxs[:, t-1] + linsystems[k].B*0.0 + linsystems[k].b

      if (-2.5 < temp[1] < 1.5) && (0.0 < temp[2] < 1000) # some of the linearisations are unbounded
        linxs[:, t] = temp
      else
        linxs[:, t] = linxs[:, t-1]
      end
  end

  subplot(nDD, 1, dd)
  x1, = plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
  x2, = plot(linxs[1,:][:], linxs[2,:][:], "r--", linewidth=3)
  plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
  plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
  plot(linxs[1,1], linxs[2,1], "ro", markersize=10, markeredgewidth = 4)
  plot(linxs[1,end], linxs[2,end], "rx", markersize=10, markeredgewidth = 4)
  ylabel("Temperature [K]")

  ## Comment out as necessary!
  # ss1 = readcsv("ss1.csv")
  # ss2 = readcsv("ss2.csv")
  ss3 = readcsv("ss3.csv")
  # x3, = plot(ss1[1], ss1[2], "gx", markersize=10, markeredgewidth = 4)
  # x3, = plot(ss2[1], ss2[2], "gx", markersize=10, markeredgewidth = 4)
  x3, = plot(ss3[1], ss3[2], "gx", markersize=10, markeredgewidth = 4)

end

legend([x1, x2, x3],["Nonlinear Model","Linear Model","Operating Point"], loc="best")
xlabel(L"Concentration [kmol.m$^{-3}$]")
rc("font",size=22)
