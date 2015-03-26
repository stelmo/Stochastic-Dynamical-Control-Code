# Linearisation Procedure
using PyPlot
import Reactor_functions

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
nX = 4 # rows
nY = 4 # cols
total_ops = nX*nY + 3 # nX*nY +3
npoints = 10
total_ops = npoints + 3 # nX*nY +3

xspace = [0.0, 1.0]
yspace = [250, 550]

# linsystems = Reactor_functions.getLinearSystems(nX, nY, xspace, yspace, h, cstr)
linsystems = Reactor_functions.getLinearSystems_randomly(npoints, xspace, yspace, h, cstr)
rc("font", family="serif", size=24)
figure(1)
x1 = 0.0
x2 = 0.0
for k=1:(total_ops)
  if k == npoints + 1
    initial_states = linsystems[k].op + abs(randn(2)).*[0.15, 20.0]
  elseif k == npoints + 2
    initial_states = linsystems[k].op + randn(2).*[0.1, 10.0]
  elseif k == npoints + 3
    initial_states = linsystems[k].op - abs(randn(2)).*[0.2, 20.0]
  else
    initial_states = linsystems[k].op + randn(2).*[0.03, 4.0]
  end

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

  x1, = plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
  x2, = plot(linxs[1,:][:], linxs[2,:][:], "r--", linewidth=3)
  plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
  plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
  plot(linxs[1,1], linxs[2,1], "ro", markersize=10, markeredgewidth = 4)
  plot(linxs[1,end], linxs[2,end], "rx", markersize=10, markeredgewidth = 4)
end
# Also plot the steady state points
ss1 = readcsv("ss1.csv")
ss2 = readcsv("ss2.csv")
ss3 = readcsv("ss3.csv")
x3, = plot(ss1[1], ss1[2], "gx", markersize=10, markeredgewidth = 4)
plot(ss2[1], ss2[2], "gx", markersize=10, markeredgewidth = 4)
plot(ss3[1], ss3[2], "gx", markersize=10, markeredgewidth = 4)
legend([x1, x2, x3],["Nonlinear Model","Linear Model","Operating Point"], loc="best")

xlim([-0.1, 1.2])
ylim([250,850])
ylabel("Temperature [K]")
xlabel(L"Concentration [kmol.m$^{-3}$]")
