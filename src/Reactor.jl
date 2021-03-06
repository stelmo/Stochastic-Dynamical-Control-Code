#  CSTR Model
# One measured variable, one manipulated variable, two states.
# Taken from the paper "Stabilization of nonlinear systems with state and control constraints using
# Lyapunov-based predictive control" by Mhaskar, P and El-Farra, N.H. and Christofides, P.D. See
# the folder: literature/control on dropbox
module Reactor

using NLsolve

type reactor
  # Supplies the parameters of the reactor. They may change i.e. drift...
  V :: Float64
  R :: Float64
  CA0:: Float64
  TA0:: Float64
  dH :: Float64
  k0 :: Float64
  E :: Float64
  Cp :: Float64
  rho :: Float64
  F :: Float64
end

type LinearReactor
  op:: Array{Float64, 1} # linearisation points
  A :: Array{Float64, 2}
  B :: Array{Float64, 1}
  b :: Array{Float64, 1} # deviation variables
end

function run_reactor(xprev::Array{Float64, 1}, u::Float64, h::Float64, model::reactor)
  # Use Runga-Kutta method to solve for the next time step using the full
  # nonlinear model.
  k1 :: Array{Float64, 1} = reactor_ode(xprev, u, model)
  k2 :: Array{Float64, 1} = reactor_ode(xprev + 0.5*h.*k1, u, model)
  k3 :: Array{Float64, 1} = reactor_ode(xprev + 0.5*h.*k2, u, model)
  k4 :: Array{Float64, 1} = reactor_ode(xprev + h.*k3, u, model)
  xnow :: Array{Float64, 1} = xprev + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
  return xnow
end

function reactor_ode(xprev::Array{Float64, 1}, u::Float64, model::reactor)
  # Evaluate the ODE functions describing the reactor.
  xnow :: Array{Float64, 1} = zeros(2)
  xnow[1] = (model.F/model.V) * (model.CA0 - xprev[1]) - model.k0*exp(-model.E/(model.R*xprev[2]))*xprev[1]
  xnow[2] = (model.F/model.V) * (model.TA0 - xprev[2]) - (model.dH/(model.rho*model.Cp))*model.k0*exp(-model.E/(model.R*xprev[2]))*xprev[1] + u/(model.rho*model.Cp*model.V)
  return xnow
end

function reactor_func!(xprev::Array{Float64, 1}, u::Float64, model::reactor, xnow::Array{Float64, 1})
  # Evaluate the ODE functions describing the reactor. In the format required by NLsolve!
  # xnow :: Array{Float64, 1} = zeros(2)
  xnow[1] = (model.F/model.V) * (model.CA0 - xprev[1]) - model.k0*exp(-model.E/(model.R*xprev[2]))*xprev[1]
  xnow[2] = (model.F/model.V) * (model.TA0 - xprev[2]) - (model.dH/(model.rho*model.Cp))*model.k0*exp(-model.E/(model.R*xprev[2]))*xprev[1] + u/(model.rho*model.Cp*model.V)
end

function jacobian(x::Array{Float64, 1}, model::reactor)
  # Returns the Jacobian evaluated at x
  J11 = -model.F/model.V-model.k0*exp(-model.E/(model.R*x[2]))
  J12 = -x[1]*model.k0*exp(-model.E/(model.R*x[2]))*(model.E/(model.R*x[2]^2))
  J21 = -model.dH/(model.rho*model.Cp)*model.k0*exp(-model.E/(model.R*x[2]))
  J22 = -(model.F/model.V + model.dH/(model.rho*model.Cp)*model.k0*exp(-model.E/(model.R*x[2]))*(model.E/(model.R*x[2]^2))*x[1])
  return [J11 J12;J21 J22]
end

function QG(T::Float64, model::reactor)
  # Return the evaluated heat generation term.
  ca = model.F/model.V*model.CA0/(model.F/model.V + model.k0*exp(-model.E/(model.R*T)))
  qg = -model.dH/(model.rho*model.Cp)*model.k0*exp(-model.E/(model.R*T))*ca
  return qg
end

function CA(T::Float64, model::reactor)
  ca = model.F/model.V*model.CA0/(model.F/model.V + model.k0*exp(-model.E/(model.R*T)))
  return ca
end

function QR(T::Float64, Q::Float64, model::reactor)
  # Return the evaluated heat removal term.
  qr = - model.F/model.V*(model.TA0 - T) - Q/(model.rho * model.V * model.Cp)
  return qr
end

function linearise(linpoint::Array{Float64, 1}, h::Float64, model::reactor)
  # Returns the linearised coefficients of the Runge Kutta method
  # given a linearisation point.
  # To solve use x(k+1) =  Ax(k) + Bu(k)

  B11 = 0.0
  B21 = 1.0/(model.rho*model.V*model.Cp)
  B = [B11;B21]
  A = jacobian(linpoint, model)
  F0 = reactor_ode(linpoint, 0.0, model) # u = 0 because Bs account for the control term
  D = A*linpoint
  # now we have x' = Ax + Bu + F0 - D = F(x)
  # now write ito deviation variables!
  newb = A\(D-F0)
  # now we have xp' = Axp + Bu where x = xp + newb

  n, = size(A)
  # Uses the bilinear transform aka the Tustin transform... google it...
  newA = (eye(n) + 0.5*A*h)*inv(eye(n) - 0.5*A*h)
  newB = inv(A)*(newA-eye(n))*B

  return newA, newB, newb
end

function discretise(nX, nY, xspace, yspace)
  # Discrete the state space into nX*nY regions.
  dx = (xspace[2] - xspace[1])/nX
  dy = (yspace[2] - yspace[1])/nY

  operatingpoints = zeros(2,nX*nY + 3) # add the three nominal points

  k = 1 #counter
  for x=1:nX
    xnow = dx*(x-1) + xspace[1] + dx*0.5
    for y=1:nY
      ynow = dy*(y-1) + yspace[1] + dy*0.5
      operatingpoints[:, k] = [xnow, ynow]
      k += 1
    end
  end
  operatingpoints[:, k] = [0.009718824131074055, 508.0562351737852]
  operatingpoints[:, k+1] = [0.48934869384879404, 412.1302612302412]
  operatingpoints[:, k+2] = [0.9996453064079288, 310.07093871841454]
  return operatingpoints
end

function discretise_randomly(npoints, xspace, yspace)
  # Perform the same action as discretise() except pick points to
  # discretise around at random.
  operatingpoints = zeros(2, npoints+3)
  if npoints==0
    k = 0
  else
    k=1
  end
  for k=1:npoints
    nx = rand()
    ny = rand()
    xnow = xspace[1] + nx*(xspace[2] - xspace[1])
    ynow = yspace[1] + ny*(yspace[2] - yspace[1])
    operatingpoints[:, k] = [xnow, ynow]
  end

  operatingpoints[:, k+1] = [0.009718824131074055, 508.0562351737852]
  operatingpoints[:, k+2] = [0.48934869384879404, 412.1302612302412]
  operatingpoints[:, k+3] = [0.9996453064079288, 310.07093871841454]
  return operatingpoints
end

function getLinearSystems(nX, nY, xspace, yspace, h, model::reactor)
  # Returns an array of linearised systems

  N = nX*nY + 3 # add the three nominal operating points
  linsystems = Array(LinearReactor, N)
  ops = discretise(nX, nY, xspace, yspace)
  for k=1:N
    op = ops[:, k]
    A, B, b = linearise(op, h, model)
    linsystems[k] = LinearReactor(op, A, B, b)
  end

  return linsystems
end

function getLinearSystems_randomly(npoints, xspace, yspace, h, model::reactor)
  # Returns an array of linearised systems

  N = npoints + 3 # add the three nominal operating points
  linsystems = Array(LinearReactor, N)
  ops = discretise_randomly(npoints, xspace, yspace)
  for k=1:N
    op = ops[:, k]
    A, B, b = linearise(op, h, model)
    linsystems[k] = LinearReactor(op, A, B, b)
  end

  return linsystems
end

function getNominalLinearSystems(h, model::reactor)
  # Returns an array of linearised systems

  linsystems = Array(LinearReactor, 3)

  #Get the steady state points
  xguess1 = [0.073, 493.0]
  xguess2 = [0.21, 467.0]
  xguess3 = [0.999, 310.0]
  f!(x, xvec) = Reactor.reactor_func!(x, 0.0, model, xvec)

  xx1res = nlsolve(f!, xguess1)
  xx2res = nlsolve(f!, xguess2)
  xx3res = nlsolve(f!, xguess3)

  ops = zeros(length(xx1res.zero), 3)
  ops[:, 1] = xx1res.zero
  ops[:, 2] = xx2res.zero
  ops[:, 3] = xx3res.zero

  for k=1:3
    op = ops[:, k]
    A, B, b = linearise(op, h, model)
    linsystems[k] = LinearReactor(op, A, B, b)
  end

  return linsystems
end


end # module
