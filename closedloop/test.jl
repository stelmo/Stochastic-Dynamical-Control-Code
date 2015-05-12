include("../params.jl")

using JuMP
using Ipopt

kfmean = [0.5, 450]
limu = 15000.0
limustep = 1000.0
horizon = 10
ysp = 0.48

m = Model(solver=IpoptSolver(print_level=0)) #

@defVar(m, x[1:2, 1:horizon])

if limu == 0.0
  @defVar(m, u[1:horizon-1])
else
  @defVar(m, -limu <= u[1:horizon-1] <= limu)
end

@addConstraint(m, x[1, 1] == kfmean[1])
@addConstraint(m, x[2, 1] == kfmean[2])


@defNLExpr(x1, (cstr_model.F/cstr_model.V) * (cstr_model.CA0 - kfmean[1]) - cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*kfmean[2]))*kfmean[1])
@defNLExpr(x2, (cstr_model.F/cstr_model.V) * (cstr_model.TA0 - kfmean[2]) - (cstr_model.dH/(cstr_model.rho*cstr_model.Cp))*cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*kfmean[2]))*kfmean[1] + u[1]/(cstr_model.rho*cstr_model.Cp*cstr_model.V))

@defVar(m, k1[1:2])
@defVar(m, k2[1:2])
@defVar(m, k3[1:2])
@defVar(m, k4[1:2])

@addNLConstraint(m, k1[1] == x1)
@addNLConstraint(m, k1[2] == x2)


@addConstraint(m, x[1, 2] == kfmean[1] + h*k1[1])
@addConstraint(m, x[2, 2] == kfmean[2] + h*k1[2])

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

@setObjective(m, Min, sum{QQ[1]*x[1, i]^2 - 2.0*ysp*QQ[1]*x[1, i] + RR*u[i]^2, i=1:horizon-1} + QQ[1]*x[1, horizon]^2 - 2.0*QQ[1]*ysp*x[1, horizon])

status = solve(m)
