module MPC

warn("MPC is hardcoded for the CSTR!")

using JuMP
using Ipopt # the back up optimisation routine
using Mosek # can't handle the nonconvex constraint

function mpc_mean(adjmean, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp, limu, limstepu, revconstr, d=zeros(2))
  # return the MPC control input using a linear system

  # m = Model(solver=IpoptSolver(print_level=0)) # chooses optimiser by itself
  m = Model(solver=MosekSolver(LOG=0)) # chooses optimiser by itself

  @defVar(m, x[1:2, 1:horizon])

  if limu == 0.0
    @defVar(m, u[1:horizon-1])
  else
    @defVar(m, -limu <= u[1:horizon-1] <= limu)
  end

  @addConstraint(m, x[1, 1] == adjmean[1])
  @addConstraint(m, x[2, 1] == adjmean[2])
  @addConstraint(m, x[1, 2] == A[1,1]*adjmean[1] + A[1,2]*adjmean[2] + B[1]*u[1] + d[1])
  @addConstraint(m, x[2, 2] == A[2,1]*adjmean[1] + A[2,2]*adjmean[2] + B[2]*u[1] + d[2])

  for k=3:horizon
    @addConstraint(m, x[1, k] == A[1,1]*x[1, k-1] + A[1,2]*x[2, k-1] + B[1]*u[k-1] + d[1])
    @addConstraint(m, x[2, k] == A[2,1]*x[1, k-1] + A[2,2]*x[2, k-1] + B[2]*u[k-1] + d[2])
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

  for k=2:horizon-1
    @addConstraint(m, u[k]-u[k-1] <= limstepu)
    @addConstraint(m, u[k]-u[k-1] >= -limstepu)
  end

  @setObjective(m, Min, sum{QQ[1]*x[1, i]^2 - 2.0*ysp*QQ[1]*x[1, i] + RR*u[i]^2 - 2.0*usp*RR*u[i], i=1:horizon-1} + QQ[1]*x[1, horizon]^2 - 2.0*QQ[1]*ysp*x[1, horizon])

  status = solve(m)

  if status != :Optimal
    warn("Mosek did not converge. Attempting to use Ipopt...")
    unow = mpc_mean_i(adjmean, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp, limu, limstepu, revconstr, d)
    return unow
  end

  return getValue(u[1]) # get the controller input
end

function mpc_mean_i(adjmean, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp, limu, limstepu, revconstr, d=zeros(2))
  # Back up optimiser.

  m = Model(solver=IpoptSolver(print_level=0)) # chooses optimiser by itself

  @defVar(m, x[1:2, 1:horizon])

  if limu == 0.0
    @defVar(m, u[1:horizon-1])
  else
    @defVar(m, -limu <= u[1:horizon-1] <= limu)
  end

  @addConstraint(m, x[1, 1] == adjmean[1])
  @addConstraint(m, x[2, 1] == adjmean[2])
  @addConstraint(m, x[1, 2] == A[1,1]*adjmean[1] + A[1,2]*adjmean[2] + B[1]*u[1] + d[1])
  @addConstraint(m, x[2, 2] == A[2,1]*adjmean[1] + A[2,2]*adjmean[2] + B[2]*u[1] + d[2])

  for k=3:horizon
    @addConstraint(m, x[1, k] == A[1,1]*x[1, k-1] + A[1,2]*x[2, k-1] + B[1]*u[k-1] + d[1])
    @addConstraint(m, x[2, k] == A[2,1]*x[1, k-1] + A[2,2]*x[2, k-1] + B[2]*u[k-1] + d[2])
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

  for k=2:horizon-1
    @addConstraint(m, u[k]-u[k-1] <= limstepu)
    @addConstraint(m, u[k]-u[k-1] >= -limstepu)
  end

  @setObjective(m, Min, sum{QQ[1]*x[1, i]^2 - 2.0*ysp*QQ[1]*x[1, i] + RR*u[i]^2 - 2.0*usp*RR*u[i], i=1:horizon-1} + QQ[1]*x[1, horizon]^2 - 2.0*QQ[1]*ysp*x[1, horizon])

  status = solve(m)

  if status == :Optimal
    info("Ipopt converged!")
  end

  return getValue(u[1]) # get the controller input
end

function mpc_var(adjmean, fcovar, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp, limu, limstepu, revconstr, swapcon, Q, sigma, growvar, d=zeros(2))
  # return the MPC control input using a linear system

  # m = Model(solver=IpoptSolver(print_level=0)) # chooses optimiser by itself
  m = Model(solver=MosekSolver(LOG=0)) # chooses optimiser by itself

  @defVar(m, x[1:2, 1:horizon])

  if limu == 0.0
    @defVar(m, u[1:horizon-1])
  else
    @defVar(m, -limu <= u[1:horizon-1] <= limu)
  end


  @addConstraint(m, x[1, 1] == adjmean[1])
  @addConstraint(m, x[2, 1] == adjmean[2])
  @addConstraint(m, x[1, 2] == A[1,1]*adjmean[1] + A[1,2]*adjmean[2] + B[1]*u[1] + d[1])
  @addConstraint(m, x[2, 2] == A[2,1]*adjmean[1] + A[2,2]*adjmean[2] + B[2]*u[1] + d[2])

  for k=3:horizon
    @addConstraint(m, x[1, k] == A[1,1]*x[1, k-1] + A[1,2]*x[2, k-1] + B[1]*u[k-1] + d[1])
    @addConstraint(m, x[2, k] == A[2,1]*x[1, k-1] + A[2,2]*x[2, k-1] + B[2]*u[k-1] + d[2])
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

  # # add distribution constraints (Ipopt doens't like it when its a quadratic constraint because nonconvex)
  # sigma = 2.2788 # one sigma 68 % confidence
  # sigma = 4.605 # 90 % confidence
  # sigma = 9.21 # 99 % confidence

  if growvar
    predvar = Q + A*fcovar*transpose(A)
    for k=2:horizon # don't do anything about k=1
      rsquared = [aline, bline]'*predvar*[aline, bline]
      r = sqrt(sigma*rsquared[1])
      @addConstraint(m, swapcon*(cline + aline*(x[1, k] + b[1]) + bline*(x[2, k] + b[2])) >= r)
      predvar = Q + A*predvar*transpose(A)
    end
  else
    predvar = Q + A*fcovar*transpose(A)
    for k=2:horizon # don't do anything about k=1
      rsquared = [aline, bline]'*predvar*[aline, bline]
      r = sqrt(sigma*rsquared[1])
      @addConstraint(m, swapcon*(cline + aline*(x[1, k] + b[1]) + bline*(x[2, k] + b[2])) >= r)
    end
  end

  for k=2:horizon-1
    @addConstraint(m, u[k]-u[k-1] <= limstepu)
    @addConstraint(m, u[k]-u[k-1] >= -limstepu)
  end

  @setObjective(m, Min, sum{QQ[1]*x[1, i]^2 - 2.0*ysp*QQ[1]*x[1, i] + RR*u[i]^2 - 2.0*usp*RR*u[i], i=1:horizon-1} + QQ[1]*x[1, horizon]^2 - 2.0*QQ[1]*ysp*x[1, horizon])

  status = solve(m)
  if status != :Optimal
    warn("Mosek did not converge. Attempting to use Ipopt...")
    unow = mpc_var_i(adjmean, fcovar, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp, limu, limstepu, revconstr, swapcon, Q, sigma, growvar, d)
    return unow
  end

  return getValue(u[1]) # get the controller input
end

function mpc_var_i(adjmean, fcovar, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp, limu, limstepu, revconstr, swapcon, Q, sigma, growvar, d=zeros(2))
  # return the MPC control input using a linear system

  m = Model(solver=IpoptSolver(print_level=0)) # chooses optimiser by itself

  @defVar(m, x[1:2, 1:horizon])

  if limu == 0.0
    @defVar(m, u[1:horizon-1])
  else
    @defVar(m, -limu <= u[1:horizon-1] <= limu)
  end


  @addConstraint(m, x[1, 1] == adjmean[1])
  @addConstraint(m, x[2, 1] == adjmean[2])
  @addConstraint(m, x[1, 2] == A[1,1]*adjmean[1] + A[1,2]*adjmean[2] + B[1]*u[1] + d[1])
  @addConstraint(m, x[2, 2] == A[2,1]*adjmean[1] + A[2,2]*adjmean[2] + B[2]*u[1] + d[2])

  for k=3:horizon
    @addConstraint(m, x[1, k] == A[1,1]*x[1, k-1] + A[1,2]*x[2, k-1] + B[1]*u[k-1] + d[1])
    @addConstraint(m, x[2, k] == A[2,1]*x[1, k-1] + A[2,2]*x[2, k-1] + B[2]*u[k-1] + d[2])
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

  for k=2:horizon-1
    @addConstraint(m, u[k]-u[k-1] <= limstepu)
    @addConstraint(m, u[k]-u[k-1] >= -limstepu)
  end

  # add distribution constraints
  # sigma = 2.2788 # one sigma 68 % confidence
  # sigma = 4.605 # 90 % confidence
  # sigma = 9.21 # 99 % confidence
  if growvar
    predvar = Q + A*fcovar*transpose(A)
    c = aline*b[1] + bline*b[2] + cline
    for k=2:horizon # don't do anything about k=1
      rsquared = [aline, bline]'*predvar*[aline, bline]
      r = rsquared[1]
      @addNLConstraint(m, aline^2*x[1, k]^2 + bline^2*x[2, k]^2 + c^2 + 2.0*aline*c*x[1, k] + 2.0*aline*bline*x[1, k]*x[2, k] + 2.0*bline*c*x[2, k] >= sigma*r)
      predvar = Q + A*predvar*transpose(A)
    end
  else
    predvar = Q + A*fcovar*transpose(A)
    c = aline*b[1] + bline*b[2] + cline
    for k=2:horizon # don't do anything about k=1
      rsquared = [aline, bline]'*predvar*[aline, bline]
      r = rsquared[1]
      @addNLConstraint(m, aline^2*x[1, k]^2 + bline^2*x[2, k]^2 + c^2 + 2.0*aline*c*x[1, k] + 2.0*aline*bline*x[1, k]*x[2, k] + 2.0*bline*c*x[2, k] >= sigma*r)
    end
  end

  @setObjective(m, Min, sum{QQ[1]*x[1, i]^2 - 2.0*ysp*QQ[1]*x[1, i] + RR*u[i]^2 - 2.0*usp*RR*u[i], i=1:horizon-1} + QQ[1]*x[1, horizon]^2 - 2.0*QQ[1]*ysp*x[1, horizon])

  status = solve(m)

  if status == :Optimal
    info("Ipopt has converged!")
  end

  return getValue(u[1]) # get the controller input
end

function mpc_lqr(adjmean, horizon, A, B, b, QQ, RR, ysp, usp, d=zeros(2))
  # return the MPC control input using a linear system

  # m = Model(solver=IpoptSolver(print_level=0)) # chooses optimiser by itself
  m = Model(solver=MosekSolver(LOG=0)) # chooses optimiser by itself

  @defVar(m, x[1:2, 1:horizon])
  @defVar(m, u[1:horizon-1])

  @addConstraint(m, x[1, 1] == adjmean[1])
  @addConstraint(m, x[2, 1] == adjmean[2])
  @addConstraint(m, x[1, 2] == A[1,1]*adjmean[1] + A[1,2]*adjmean[2] + B[1]*u[1] + d[1])
  @addConstraint(m, x[2, 2] == A[2,1]*adjmean[1] + A[2,2]*adjmean[2] + B[2]*u[1] + d[2])

  for k=3:horizon
    @addConstraint(m, x[1, k] == A[1,1]*x[1, k-1] + A[1,2]*x[2, k-1] + B[1]*u[k-1] + d[1])
    @addConstraint(m, x[2, k] == A[2,1]*x[1, k-1] + A[2,2]*x[2, k-1] + B[2]*u[k-1] + d[2])
  end


  @setObjective(m, Min, sum{QQ[1]*x[1, i]^2 - 2.0*ysp*QQ[1]*x[1, i] + RR*u[i]^2 - 2.0*usp*RR*u[i], i=1:horizon-1} + QQ[1]*x[1, horizon]^2 - 2.0*QQ[1]*ysp*x[1, horizon])

  status = solve(m)

  if status != :Optimal
    warn("Mosek did not converge. Attempting to use Ipopt...")
    unow = mpc_lqr_i(adjmean, horizon, A, B, b, QQ, RR, ysp, usp, d)
    return unow
  end

  return getValue(u[1]) # get the controller input
end

function mpc_lqr_i(adjmean, horizon, A, B, b, QQ, RR, ysp, usp, d=zeros(2))
  # Back up optimiser.

  m = Model(solver=IpoptSolver(print_level=0)) # chooses optimiser by itself

  @defVar(m, x[1:2, 1:horizon])
  @defVar(m, u[1:horizon-1])

  @addConstraint(m, x[1, 1] == adjmean[1])
  @addConstraint(m, x[2, 1] == adjmean[2])
  @addConstraint(m, x[1, 2] == A[1,1]*adjmean[1] + A[1,2]*adjmean[2] + B[1]*u[1] + d[1])
  @addConstraint(m, x[2, 2] == A[2,1]*adjmean[1] + A[2,2]*adjmean[2] + B[2]*u[1] + d[2])

  for k=3:horizon
    @addConstraint(m, x[1, k] == A[1,1]*x[1, k-1] + A[1,2]*x[2, k-1] + B[1]*u[k-1] + d[1])
    @addConstraint(m, x[2, k] == A[2,1]*x[1, k-1] + A[2,2]*x[2, k-1] + B[2]*u[k-1] + d[2])
  end


  @setObjective(m, Min, sum{QQ[1]*x[1, i]^2 - 2.0*ysp*QQ[1]*x[1, i] + RR*u[i]^2 - 2.0*usp*RR*u[i], i=1:horizon-1} + QQ[1]*x[1, horizon]^2 - 2.0*QQ[1]*ysp*x[1, horizon])

  status = solve(m)

  if status == :Optimal
    info("Ipopt converged!")
  end

  return getValue(u[1]) # get the controller input
end

end
