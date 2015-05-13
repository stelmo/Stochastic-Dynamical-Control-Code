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
nH = 5 #horizon

m = Model(solver=IpoptSolver(print_level=0)) #

@defVar(m, x1[1:nP, 1:nH]) # concen
@defVar(m, x2[1:nP, 1:nH]) # temp

if limu == 0.0
  @defVar(m, u[1:nH-1])
else
  @defVar(m, -limu <= u[1:nH-1] <= limu)
end

for np=1:nP # initial state for all particles
  @addConstraint(m, x1[np, 1] == particles.x[1, np])
  @addConstraint(m, x2[np, 1] == particles.x[2, np])
end

#
for np=1:nP
  @addNLConstraint(m, x1[np, 2] == x1[np, 1] + h*((cstr_model.F/cstr_model.V) * (cstr_model.CA0 - x1[np, 2]) - cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*x2[np, 2]))*x1[np, 2]))
#
# @addNLConstraint(m, x2[np, nh+1] == x2[np, nh] + (cstr_model.F/cstr_model.V) * (cstr_model.TA0 - x2[np, nh]) - (cstr_model.dH/(cstr_model.rho*cstr_model.Cp))*cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*x2[np, nh]))*x1[np, nh] + u[nh]/(cstr_model.rho*cstr_model.Cp*cstr_model.V))
end

# @defNLExpr(k11expr[np=1:nP, nh=1:nH-1], (cstr_model.F/cstr_model.V) * (cstr_model.CA0 - x1[np, nh]) - cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*x2[np, nh]))*x1[np, nh])
#
# @defNLExpr(k12expr[np=1:nP, nh=1:nH-1], (cstr_model.F/cstr_model.V) * (cstr_model.TA0 - x2[np, nh]) - (cstr_model.dH/(cstr_model.rho*cstr_model.Cp))*cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*x2[np, nh]))*x1[np, nh] + u[nh]/(cstr_model.rho*cstr_model.Cp*cstr_model.V))

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

#     @addNLConstraint(m, x1[np, nh+1] == x1[np, nh] + (cstr_model.F/cstr_model.V) * (cstr_model.CA0 - x1[np, nh]) - cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*x2[np, nh]))*x1[np, nh])
#     @addNLConstraint(m, x2[np, nh+1] == x2[np, nh] + (cstr_model.F/cstr_model.V) * (cstr_model.TA0 - x2[np, nh]) - (cstr_model.dH/(cstr_model.rho*cstr_model.Cp))*cstr_model.k0*exp(-cstr_model.E/(cstr_model.R*x2[np, nh]))*x1[np, nh] + u[nh]/(cstr_model.rho*cstr_model.Cp*cstr_model.V))
#
#   end
# end
#

@setObjective(m, Min, sum{u[i]^2, i=1:nH-1})


#
status = solve(m)
