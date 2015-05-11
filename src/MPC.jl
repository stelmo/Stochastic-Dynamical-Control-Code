module MPC

warn("MPC is hardcoded for the CSTR!")

using JuMP
using Ipopt

function mpc_mean(kfmean, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, limu, revconstr)
  # return the MPC control input using a linear system

  m = Model(solver=IpoptSolver(print_level=0)) # chooses optimiser by itself

  @defVar(m, x[1:2, 1:horizon])

  if limu == 0.0
    @defVar(m, u[1:horizon-1])
  else
    @defVar(m, -limu <= u[1:horizon-1] <= limu)
  end

  @addConstraint(m, x[1, 1] == kfmean[1])
  @addConstraint(m, x[2, 1] == kfmean[2])
  @addConstraint(m, x[1, 2] == A[1,1]*kfmean[1] + A[1,2]*kfmean[2] + B[1]*u[1])
  @addConstraint(m, x[2, 2] == A[2,1]*kfmean[1] + A[2,2]*kfmean[2] + B[2]*u[1])

  for k=3:horizon
    @addConstraint(m, x[1, k] == A[1,1]*x[1, k-1] + A[1,2]*x[2, k-1] + B[1]*u[k-1])
    @addConstraint(m, x[2, k] == A[2,1]*x[1, k-1] + A[2,2]*x[2, k-1] + B[2]*u[k-1])
  end

  # # add state constraints
  if revconstr
    for k=2:horizon # can't do anything about k=1
      @addConstraint(m, aline*(x[1, k] + b[1]) + bline*(x[2, k] + b[2]) <= -1.0*cline)
    end
  else
    for k=2:horizon # can't do anything about k=1
      @addConstraint(m, aline*(x[1, k] + b[1]) + bline*(x[2, k] + b[2]) >= -1.0*cline)
    end
  end


  @setObjective(m, Min, sum{QQ[1]*x[1, i]^2 - 2.0*ysp*QQ[1]*x[1, i] + RR*u[i]^2, i=1:horizon-1} + QQ[1]*x[1, horizon]^2 - 2.0*QQ[1]*ysp*x[1, horizon])

  status = solve(m)

  return getValue(u[1]) # get the controller input
end

function mpc_var(kfmean, kfcovar, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, limu)
  # return the MPC control input using a linear system

  m = Model(solver=IpoptSolver(print_level=0)) # chooses optimiser by itself

  @defVar(m, x[1:2, 1:horizon])

  if limu
    @defVar(m, -15000.0 <= u[1:horizon-1] <= 15000.0)
  else
    @defVar(m, u[1:horizon-1])
  end

  @addConstraint(m, x[1, 1] == kfmean[1])
  @addConstraint(m, x[2, 1] == kfmean[2])
  @addConstraint(m, x[1, 2] == A[1,1]*kfmean[1] + A[1,2]*kfmean[2] + B[1]*u[1])
  @addConstraint(m, x[2, 2] == A[2,1]*kfmean[1] + A[2,2]*kfmean[2] + B[2]*u[1])

  for k=3:horizon
    @addConstraint(m, x[1, k] == A[1,1]*x[1, k-1] + A[1,2]*x[2, k-1] + B[1]*u[k-1])
    @addConstraint(m, x[2, k] == A[2,1]*x[1, k-1] + A[2,2]*x[2, k-1] + B[2]*u[k-1])
  end

  # # add state constraints
  for k=2:horizon # can't do anything about k=1
    @addConstraint(m, aline*(x[1, k] + b[1]) + bline*(x[2, k] + b[2]) >= -1.0*cline)
  end

  # add distribution constraints
  sigma = 4.605 # this should match the chi square value in Ellipse
  for k=2:horizon # don't do anything about k=1
    @addNLConstraint(m, (cline + aline*(x[1, k] + b[1]) + bline*(x[2, k] + b[2]))^2/(aline^2*kfcovar[1,1] + bline^2*kfcovar[2,2] + aline*bline*kfcovar[1,2] + aline*bline*kfcovar[2,1]) >= sigma)
  end

  @setObjective(m, Min, sum{QQ[1]*x[1, i]^2 - 2.0*ysp*QQ[1]*x[1, i] + RR*u[i]^2, i=1:horizon-1} + QQ[1]*x[1, horizon]^2 - 2.0*QQ[1]*ysp*x[1, horizon])

  status = solve(m)

  return getValue(u[1]) # get the controller input
end


end
