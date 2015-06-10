# Check the stability of the different linearisation points over time
# Note: small time spans required here

include("./params.jl") # load all the parameters and modules

# Divide state space into sectors: n by m
nX = 20 # rows
nY = 20 # cols
total_ops = nX*nY # ignore the nominal ss points

xspace = [0.0, 1.0]
yspace = [250, 650]

linsystems = Reactor.getLinearSystems(nX, nY, xspace, yspace, h, cstr_model)

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
    linxs[:,1] = initial_states - linsystems[k].b
    # Loop through the rest of time
    flag = false
    for t=2:N
        temp1 = Reactor.run_reactor(xs[:, t-1], 0.0, h, cstr_model) # actual plant
        temp2 = linsystems[k].A*linxs[:, t-1] + linsystems[k].B*0.0
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
    linxs = linxs .+ linsystems[k].b
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
