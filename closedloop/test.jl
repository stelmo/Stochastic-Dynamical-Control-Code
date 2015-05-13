include("../params.jl")

using JuMP
using Ipopt

init_state = [0.5; 400] # initial state

f(x, u, w) = Reactor.run_reactor(x, u, h, cstr_model) + w
g(x) = C2*x # state observation

cstr_pf = PF.Model(f,g)

# Initialise the PF
nP = 1 # number of particles.
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
nH = 5 #horizon

m = Model(solver=IpoptSolver(print_level=0)) #

@defVar(m, x1[1:nP, 1:nH]) # concen
@defVar(m, x2[1:nP, 1:nH]) # temp

if limu == 0.0
  @defVar(m, u[1:nH-1])
else
  @defVar(m, -limu <= u[1:nH-1] <= limu)
end
#
# @defVar(m, k11[1:nP, 1:nH-1])
# @defVar(m, k12[1:nP, 1:nH-1])
#
# @defVar(m, k21[1:nP, 1:nH-1])
# @defVar(m, k22[1:nP, 1:nH-1])
#
# @defVar(m, k31[1:nP, 1:nH-1])
# @defVar(m, k32[1:nP, 1:nH-1])
#
# @defVar(m, k41[1:nP, 1:nH-1])
# @defVar(m, k42[1:nP, 1:nH-1])


for np=1:nP # initial state for all particles
  @addConstraint(m, x1[np, 1] == particles.x[1, np])
  @addConstraint(m, x2[np, 1] == particles.x[2, np])
end

# @defNLExpr(k11expr[np=1:nP, nh=1:nH-1], (cstr_model.F/cstr_model.V) * (cstr_model.CA0 - x1[np, nh]) - cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*x2[np, nh]))*x1[np, nh])
#
# @defNLExpr(k12expr[np=1:nP, nh=1:nH-1], (cstr_model.F/cstr_model.V) * (cstr_model.TA0 - x2[np, nh]) - (cstr_model.dH/(cstr_model.rho*cstr_model.Cp))*cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*x2[np, nh]))*x1[np, nh] + u[nh]/(cstr_model.rho*cstr_model.Cp*cstr_model.V))


## loop over horizon
for nh=1:nH-1
  for np=1:nP # over each particle

    # @addNLConstraint(m, k11expr[np, nh] == k11[np, nh])

    # @addNLConstraint(m, k12[np, nh] == k12expr[np, nh])

    # (x1[np, nh-1] + 0.5*h*k11[np, nh-1])
    # (x2[np, nh-1] + 0.5*h*k12[np, nh-1])
    #
    # @addNLConstraint(m, k21[np, nh-1] == (cstr_model.F/cstr_model.V) * (cstr_model.CA0 - (x1[np, nh-1] + 0.5*h*k11[np, nh-1])) - cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*(x2[np, nh-1] + 0.5*h*k12[np, nh-1])))*(x1[np, nh-1] + 0.5*h*k11[np, nh-1]))
    #
    # @addNLConstraint(m, k22[np, nh-1] == (cstr_model.F/cstr_model.V) * (cstr_model.TA0 - (x2[np, nh-1] + 0.5*h*k12[np, nh-1])) - (cstr_model.dH/(cstr_model.rho*cstr_model.Cp))*cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*(x2[np, nh-1] + 0.5*h*k12[np, nh-1])))*(x1[np, nh-1] + 0.5*h*k11[np, nh-1]) + u[nh-1]/(cstr_model.rho*cstr_model.Cp*cstr_model.V))
    #
    # @addNLConstraint(m, x2[np, nh] == x2[np, nh-1] + (h/6.0)*(k12[np, nh-1] + 2.0*k21[np, nh-1] + 2.0*k32[np, nh-1] + k42[np, nh-1]))
    # @addNLConstraint(m, x1[np, nh] == x1[np, nh-1] + (h/6.0)*(k11[np, nh-1] + 2.0*k21[np, nh-1] + 2.0*k31[np, nh-1] + k41[np, nh-1]))

    @addNLConstraint(m, x1[np, nh+1] == x1[np, nh] + (cstr_model.F/cstr_model.V) * (cstr_model.CA0 - x1[np, nh]) - cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*x2[np, nh]))*x1[np, nh])
    @addNLConstraint(m, x2[np, nh+1] == x2[np, nh] + (cstr_model.F/cstr_model.V) * (cstr_model.TA0 - x2[np, nh]) - (cstr_model.dH/(cstr_model.rho*cstr_model.Cp))*cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*x2[np, nh]))*x1[np, nh] + u[nh]/(cstr_model.rho*cstr_model.Cp*cstr_model.V))

  end
end

# end setting constraints over the horizon

# Reactor.run_reactor(xs1[:, t-1], us[t-1], h, cstr_model)[1]
# function run_reactor(xprev::Array{Float64, 1}, u::Float64, h::Float64, model::reactor)
#   # Use Runga-Kutta method to solve for the next time step using the full
#   # nonlinear model.
#   k1 :: Array{Float64, 1} = reactor_ode(xprev, u, model)
#   k2 :: Array{Float64, 1} = reactor_ode(xprev + 0.5*h.*k1, u, model)
#   k3 :: Array{Float64, 1} = reactor_ode(xprev + 0.5*h.*k2, u, model)
#   k4 :: Array{Float64, 1} = reactor_ode(xprev + h.*k3, u, model)
#   xnow :: Array{Float64, 1} = xprev + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
#   return xnow
# end
#
#
# xnow[1] = (model.F/model.V) * (model.CA0 - xprev[1]) - model.k0*exp(-model.E/(model.R*xprev[2]))*xprev[1]
# xnow[2] = (model.F/model.V) * (model.TA0 - xprev[2]) - (model.dH/(model.rho*model.Cp))*model.k0*exp(-model.E/(model.R*xprev[2]))*xprev[1] + u/(model.rho*model.Cp*model.V)


# @addConstraint(m, x[1, 2] == Reactor.run_reactor(kfmeans, u[1], h, cstr_model)[1])
# @addConstraint(m, x[2, 2] == Reactor.run_reactor(kfmeans, u[1], h, cstr_model)[2])

# for k=3:horizon
#   @addConstraint(m, x[1, k] == A[1,1]*x[1, k-1] + A[1,2]*x[2, k-1] + B[1]*u[k-1])
#   @addConstraint(m, x[2, k] == A[2,1]*x[1, k-1] + A[2,2]*x[2, k-1] + B[2]*u[k-1])
# end

# # add state constraints
# if revconstr
#   for k=2:horizon # can't do anything about k=1
#     @addConstraint(m, aline*(x[1, k] + b[1]) + bline*(x[2, k] + b[2]) <= -1.0*cline)
#   end
# else
#   for k=2:horizon # can't do anything about k=1
#     @addConstraint(m, aline*(x[1, k] + b[1]) + bline*(x[2, k] + b[2]) >= -1.0*cline)
#   end
# end

# for k=2:horizon-1
#   @addConstraint(m, u[k]-u[k-1] <= limstepu)
#   @addConstraint(m, u[k]-u[k-1] >= -limstepu)
# end

@setObjective(m, Min, sum{x1[1, i]^2, i=1:nH})

# @setObjective(m, Min, sum{QQ[1]*x[1, i]^2 - 2.0*ysp*QQ[1]*x[1, i] + RR*u[i]^2, i=1:horizon-1} + QQ[1]*x[1, horizon]^2 - 2.0*QQ[1]*ysp*x[1, horizon])
#
status = solve(m)
