module PMPC


warn("PMPC is hardcoded for the CSTR!")

using JuMP
using Ipopt # the back up optimisation routine
using Mosek # can't handle the nonconvex constraint
using Distributions

function mpc(adjmean, nP, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp, limu, limstepu, revconstr, randvars1, randvars2)
  # return the MPC control input using a linear system

  # m = Model(solver=IpoptSolver(print_level=0)) # chooses optimiser by itself
  m = Model(solver=MosekSolver(LOG=0)) # chooses optimiser by itself

  @defVar(m, x[1:2, 1:nP, 1:horizon])
  @defVar(m, c[1:nP, 1:horizon], Bin) # integer

  if limu == 0.0
    @defVar(m, u[1:horizon-1])
  else
    @defVar(m, -limu <= u[1:horizon-1] <= limu)
  end

  # initial
  for p=1:nP

    @addConstraint(m, x[1, p, 1] == adjmean[1, p])
    @addConstraint(m, x[2, p, 1] == adjmean[2, p])
    @addConstraint(m, x[1, p, 2] == A[1,1]*adjmean[1, p] + A[1,2]*adjmean[2, p] + B[1]*u[1] + randvars1[1, p])
    @addConstraint(m, x[2, p, 2] == A[2,1]*adjmean[1, p] + A[2,2]*adjmean[2, p] + B[2]*u[1] + randvars2[1, p])

    for k=3:horizon
      @addConstraint(m, x[1, p, k] == A[1,1]*x[1, p, k-1] + A[1,2]*x[2, p, k-1] + B[1]*u[k-1] + randvars1[k-1, p])
      @addConstraint(m, x[2, p, k] == A[2,1]*x[1, p, k-1] + A[2,2]*x[2, p, k-1] + B[2]*u[k-1] + randvars2[k-1, p])
    end

  end

  # add state constraints
  Mbig = -1e6 # some negative big number
  for p=1:nP
    for k=2:horizon # can't do anything about k=1
      @addConstraint(m, aline*(x[1, p, k] + b[1]) + bline*(x[2, p, k] + b[2]) + 1.0*cline >= Mbig*c[p, k-1])
    end
  end

  maxfrac = int(nP*0.5)
  for k=1:horizon-1
    @addConstraint(m, sum{c[i, k], i=1:nP} <= maxfrac )
  end


  for k=2:horizon-1
    @addConstraint(m, u[k]-u[k-1] <= limstepu)
    @addConstraint(m, u[k]-u[k-1] >= -limstepu)
  end

  @setObjective(m, Min, sum{QQ[1]*x[1, p, i]^2 - 2.0*ysp*QQ[1]*x[1, p, i] + RR*u[i]^2 - 2.0*usp*RR*u[i], i=1:horizon-1, p=1:nP} + sum{QQ[1]*x[1, p, horizon]^2 - 2.0*QQ[1]*ysp*x[1, p, horizon], p=1:nP})

  status = solve(m)

  # if status != :Optimal
  #   warn("Mosek did not converge. Attempting to use Ipopt...")
  #   unow = mpc_mean_i(adjmean, horizon, A, B, b, aline, bline, cline, QQ, RR, ysp, usp, limu, limstepu, revconstr)
  #   return unow
  # end

  return getValue(u[1]) # get the controller input
end

function getRandomTable(noisedist, nP, nH)
  # Get the random variable table
  rt1 = zeros(nH-1, nP)
  rt2 = zeros(nH-1, nP)

  for p=1:nP
    for k=1:nH-1
      temp = rand(noisedist)
      rt1[k, p] = temp[1]
      rt2[k, p] = temp[2]
    end
  end
  return rt1, rt2
end

end # module
