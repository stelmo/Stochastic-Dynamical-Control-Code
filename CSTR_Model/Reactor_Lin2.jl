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
tend = 5.0 # end simulation time
ts = [0.0:h:tend]

# Divide state space into sectors: n by m
nX = 20 # rows
nY = 20 # cols
total_ops = nX*nY # ignore the nominal ss points

xspace = [0.0, 1.0]
yspace = [250, 650]

linsystems = Reactor_functions.getLinearSystems(nX, nY, xspace, yspace, h, cstr)

diff = zeros(nY, nX)
xpoints = zeros(nX)
ypoints = zeros(nY)
k = 1 # counter
temp1  = zeros(2)
temp2 = zeros(2)
initial_states = zeros(2)
for x=1:nX

  for y=1:nY

    initial_states = linsystems[k].op
    N = length(ts)
    xs = zeros(2, N)
    linxs = zeros(2, N)
    xs[:,1] = initial_states
    linxs[:,1] = initial_states
    # Loop through the rest of time
    flag = false
    for t=2:N
        temp1 = Reactor_functions.run_reactor(xs[:, t-1], 0.0, h, cstr) # actual plant
        temp2 = linsystems[k].A*linxs[:, t-1] + linsystems[k].B*0.0 + linsystems[k].b
        if isnan(temp1[1]) || isnan(temp1[2])
          flag = true
          break
        else
          xs[:, t] = temp1
        end
        if isnan(temp2[1]) || isnan(temp2[2])
          flag = true
          break
        else
          linxs[:, t] = temp2
        end
    end
    if flag
      diff[y,x] = -1.0
    else
      diff[y, x] = norm((xs[:, end] - linxs[:, end])./[xs[:, end]], 1)
    end
    ypoints[y] = initial_states[2]
    k += 1
  end
  xpoints[x] = initial_states[1]
end

setmax = 1.0
for k=1:length(diff)
  if diff[k] == -1.0 || diff[k] > setmax
      diff[k] = setmax
  end
end

rc("font", family="serif", size=24)
figure(1)
contourf(xpoints, ypoints, diff, 50, cmap = "cubehelix")
xlabel(L"Concentration [kmol.m$^{-3}$]")
ylabel("Temperature [K]")
colorbar()
