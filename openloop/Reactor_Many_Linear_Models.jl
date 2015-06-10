# Compare the state space response of multiple linear models
# Note that small time spans are required here

include("./params.jl") # load all the parameters and modules

# Divide state space into sectors: n by m
nX = 4 # rows
nY = 4 # cols
total_ops = nX*nY + 3 # nX*nY +3
npoints = 10
total_ops = npoints + 3 # nX*nY +3

xspace = [0.0, 1.0]
yspace = [250, 550]

# linsystems = Reactor.getLinearSystems(nX, nY, xspace, yspace, h, cstr_model)
linsystems = Reactor.getLinearSystems_randomly(npoints, xspace, yspace, h, cstr_model)
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
  linxs[:,1] = initial_states - linsystems[k].b
  # Loop through the rest of time
  for t=2:N
      xs[:, t] = Reactor.run_reactor(xs[:, t-1], 0.0, h, cstr_model) # actual plant
      linxs[:, t] = linsystems[k].A*linxs[:, t-1] + linsystems[k].B*0.0
  end

  linxs = linxs .+ linsystems[k].b
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
