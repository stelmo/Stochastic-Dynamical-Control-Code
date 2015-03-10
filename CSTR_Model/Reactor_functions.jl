#  CSTR Model
# One measured variable, one manipulated variable, two states.
# Taken from the paper "Stabilization of nonlinear systems with state and control constraints using
# Lyapunov-based predictive control" by Mhaskar, P and El-Farra, N.H. and Christofides, P.D. See
# the folder: literature/control on dropbox
module Reactor_functions

type Reactor
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

function run_reactor(xprev::Array{Float64, 1}, u::Float64, h::Float64, model::Reactor)
  # Use Runga-Kutta method to solve for the next time step using the full
  # nonlinear model.
  k1 :: Array{Float64, 1} = reactor_ode(xprev, u, model)
  k2 :: Array{Float64, 1} = reactor_ode(xprev + 0.5*h.*k1, u, model)
  k3 :: Array{Float64, 1} = reactor_ode(xprev + 0.5*h.*k2, u, model)
  k4 :: Array{Float64, 1} = reactor_ode(xprev + h.*k3, u, model)
  xnow :: Array{Float64, 1} = xprev + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
  return xnow
end

function reactor_ode(xprev::Array{Float64, 1}, u::Float64, model::Reactor)
  # Evaluate the ODE functions describing the reactor.
  xnow :: Array{Float64, 1} = zeros(2)
  xnow[1] = (model.F/model.V) * (model.CA0 - xprev[1]) - model.k0*exp(-model.E/(model.R*xprev[2]))*xprev[1]
  xnow[2] = (model.F/model.V) * (model.TA0 - xprev[2]) - (model.dH/(model.rho*model.Cp))*model.k0*exp(-model.E/(model.R*xprev[2]))*xprev[1] + u/(model.rho*model.Cp*model.V)
  return xnow
end

function reactor_jacobian(x::Array{Float64, 1}, model::Reactor)
  # Returns the Jacobian evaluated at x
  J11 = -model.F/model.V-model.k0*exp(-model.E/(model.R*x[2]))
  J12 = -x[1]*model.k0*exp(-model.E/(model.R*x[2]))*(model.E/(model.R*x[2]^2))
  J21 = -model.dH/(model.rho*model.Cp)*model.k0*exp(-model.E/(model.R*x[2]))
  J22 = -(model.F/model.V + model.dH/(model.rho*model.Cp)*model.k0*exp(-model.E/(model.R*x[2]))*(model.E/(model.R*x[2]^2))*x[1])
end

function QG(T::Float64, model::Reactor)
  # Return the evaluated heat generation term.
  ca = model.F/model.V*model.CA0/(model.F/model.V + model.k0*exp(-model.E/(model.R*T)))
  qg = -model.dH/(model.rho*model.Cp)*model.k0*exp(-model.E/(model.R*T))*ca
end

function QR(T::Float64, Q::Float64, model::Reactor)
  # Return the evaluated heat removal term.
  qr = - model.F/model.V*(model.TA0 - T) - Q/(model.rho * model.V * model.Cp)
  return qr
end

end # module
