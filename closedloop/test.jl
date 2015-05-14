include("../params.jl")

using JuMP
using Ipopt

init_state = [0.5; 400] # initial state

f(x, u, w) = Reactor.run_reactor(x, u, h, cstr_model) + w
g(x) = C2*x # state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 5 # number of particles.
prior_dist = MvNormal(init_state, init_state_covar) # prior distribution
particles = PF.init_PF(prior_dist, nP, 2) # initialise the particles
state_noise_dist = MvNormal(Q) # state distribution
meas_noise_dist = MvNormal(R2) # measurement distribution

# Time step 1
xs[:,1] = init_state
ys2[:, 1] = C2*xs[:, 1] + rand(meas_noise_dist) # measured from actual plant
PF.init_filter!(particles, 0.0, ys2[:, 1], meas_noise_dist, cstr_pf)


limu = 15000.0
limustep = 1000.0
horizon = 10
ysp = 0.48
nH = 3 #horizon


## Start IPOPT

function eval_f(x)
  # Objective function
  totvar = 3*nH-1 #(x1, x2, u)
  ustart = 2*nH+1 # remember x = (concentration and temperature)
  return sum(x[ustart:1:end].^2)
end

function eval_grad_f(x, grad_f)
  # Gradient of the objective function
  for i=1:nH-1
    grad_f[i] = 2.0
  end
end

function eval_g(x, g)
  # Constraints
  # parameters
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

  # assume Forward Euler is good enough!
  function x1(x)
    return (F/V) * (CA0 - x[1]) - k0*exp(-E/(R*x[2]))*x[1]
  end

  function x2(x, u)
    return (F/V) * (TA0 - x[2]) - (dH/(rho*Cp))*k0*exp(-E/(R*x[2]))*x[1] + u/(rho*Cp*V)
  end

  totvar = 3*nH-1
  ustart = 2*nH+1

  g[1] = x[1] - particles.x[1,1] # needs to be zero
  g[2] = x[2] - particles.x[2,1] # needs to be zero

  ucounter = ustart
  for i=3:2:ustart-1
    g[i] = x[i] - x[i-2] - h*x1([x[i-2], x[i-1]])
    g[i+1] = x[i+1] - x[i-1] - h*x2([x[i-2], x[i-1]], x[ucounter])
    ucounter += 1
  end
end

function eval_jac_g(x, mode, rows, cols, value)
  # parameters
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

  if mode == :Structure
    rows[1] = 1; cols[1] = 1
    rows[2] = 2; cols[2] = 2

    rowcount = 3
    colown = 1
    ucount = 2*nH+1

    totdyncons = 2 + 4*(2*(nH-1))
    for i=3:8:totdyncons # in terms of horizon
      rows[i] = rowcount; cols[i] = colown
      rows[i+1] = rowcount; cols[i+1] = colown + 1
      rows[i+2] = rowcount; cols[i+2] = colown + 2
      rows[i+3] = rowcount; cols[i+3] = ucount

      rows[i+4] = rowcount + 1; cols[i+4] = colown
      rows[i+5] = rowcount + 1; cols[i+5] = colown + 1
      rows[i+6] = rowcount + 1; cols[i+6] = colown + 3
      rows[i+7] = rowcount + 1; cols[i+7] = ucount
      rowcount += 2
      ucount += 1
      colown += 2
    end

  else
    value[1] = 1.0
    value[2] = 1.0

    varval = 1
    totdyncons = 2 + 4*(2*(nH-1))
    for i=3:8:totdyncons # in terms of horizon
      value[i] = -(1.0 + h*(-F/V-k0*exp(-E/(R*x[varval+1]))))
      value[i+1] = -h*(-x[varval]*k0*exp(-E/(R*x[varval+1]))*(E/(R*x[varval+1]^2)))
      value[i+2] = 1.0
      value[i+3] = 0.0

      value[i+4] = -h*(-dH/(rho*Cp)*k0*exp(-E/(R*x[varval+1])))
      value[i+5] = -(1.0 -h*(F/V + dH/(rho*Cp)*k0*exp(-E/(R*x[varval+1]))*(E/(R*x[varval+1]^2))*x[varval]))
      value[i+6] = 1.0
      value[i+7] = 1.0/(rho*Cp*V)
      varval += 2
    end
  end
end

numvars = 3*nH-1
x_L = -1.0e5*ones(numvars)
x_U = 1.0e5*ones(numvars)

numcons = nH*2
g_L = -0.1*ones(numcons)
g_U = 0.1*ones(numcons)
totdyncons = 2 + 4*(2*(nH-1))

prob = createProblem(numvars, x_L, x_U, numcons, g_L, g_U, totdyncons, 10,
                     eval_f, eval_g, eval_grad_f, eval_jac_g)

addOption(prob, "hessian_approximation", "limited-memory")

prob.x = [particles.x[1,1], particles.x[2,1], zeros(numvars-2)]
status = solveProblem(prob)

println(Ipopt.ApplicationReturnStatus[status])
println(prob.x)
println(prob.obj_val)
