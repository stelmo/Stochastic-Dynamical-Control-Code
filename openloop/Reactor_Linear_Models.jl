# Linearisation Procedure

tend = 50
include("./params.jl") # load all the parameters and modules

xspace = [0.0, 1.0]
yspace = [250, 650]
linsystems = Reactor.getLinearSystems_randomly(0, xspace, yspace, h, cstr_model)


rc("font", family="serif", serif="Computer Modern", size=32)
rc("text", usetex=true)
figure(1)
k=3 # set which operating point to use
# also remember to change +- on line 47 and the SS points on lines 75-81
nDD = 2
x1 = 0 # legend var
x2 = 0 # legend var
x3 = 0 # legend var
for dd=1:nDD # only loop through
  initial_states = linsystems[k].op + randn(2).*[0.01;10]*dd

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

  subplot(nDD, 1, dd)
  x1, = plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
  x2, = plot(linxs[1,:][:], linxs[2,:][:], "r--", linewidth=3)
  plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
  plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
  plot(linxs[1,1], linxs[2,1], "ro", markersize=10, markeredgewidth = 4)
  plot(linxs[1,end], linxs[2,end], "rx", markersize=10, markeredgewidth = 4)
  ylabel(L"T$_R$ [K]")
  locator_params(nbins=6)


  ## Comment out as necessary!
  if k==1
    ss1 = [0.0097, 508.0562]
    x3, = plot(ss1[1], ss1[2], "gx", markersize=10, markeredgewidth = 4)
  elseif k==2
    ss2 = [0.4893, 412.1302]
    x3, = plot(ss2[1], ss2[2], "gx", markersize=10, markeredgewidth = 4)
  else
    ss3 = [0.9996, 310.0709]
    x3, = plot(ss3[1], ss3[2], "gx", markersize=10, markeredgewidth = 4)
  end
end

legend([x1, x2, x3],["Nonlinear model","Linear model","Operating point"], loc="best")
xlabel(L"C$_A$ [kmol.m$^{-3}$]")
locator_params(nbins=6)
