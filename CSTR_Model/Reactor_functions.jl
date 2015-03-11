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

function jacobian(x::Array{Float64, 1}, model::Reactor)
  # Returns the Jacobian evaluated at x
  J11 = -model.F/model.V-model.k0*exp(-model.E/(model.R*x[2]))
  J12 = -x[1]*model.k0*exp(-model.E/(model.R*x[2]))*(model.E/(model.R*x[2]^2))
  J21 = -model.dH/(model.rho*model.Cp)*model.k0*exp(-model.E/(model.R*x[2]))
  J22 = -(model.F/model.V + model.dH/(model.rho*model.Cp)*model.k0*exp(-model.E/(model.R*x[2]))*(model.E/(model.R*x[2]^2))*x[1])
  return [J11 J12;J21 J22]
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

function linearise(linpoint::Array{Float64, 1}, h::Float64, model::Reactor)
  # Returns the linearised coefficients given a linearisation point.

  B11 = 0.0
  B21 = 1.0/(model.rho*model.V*model.Cp)
  jac = jacobian(linpoint, model)
  J11 = jac[1,1]
  J12 = jac[1,2]
  J21 = jac[2,1]
  J22 = jac[2,2]
  rfunc = reactor_ode(linpoint, 0.0, model) # u = 0 because Bs account for the control term
  f0 = rfunc[1]
  g0 = rfunc[2]
  D = jac*linpoint
  D11 = D[1]
  D21 = D[2]

  A11 = ((J11^4*h^4)/24 + (J11^3*h^3)/6 + (J11^2*J12*J21*h^4)/8 + (J11^2*h^2)/2 + (J11*J12*J21*J22*h^4)/12 + (J11*J12*J21*h^3)/3 + J11*h + (J12^2*J21^2*h^4)/24 + (J12*J21*J22^2*h^4)/24 + (J12*J21*J22*h^3)/6 + (J12*J21*h^2)/2)*linpoint[1]
  A12 = ((J11^3*J12*h^4)/24 + (J11^2*J12*J22*h^4)/24 + (J11^2*J12*h^3)/6 + (J21*J11*J12^2*h^4)/12 + (J11*J12*J22^2*h^4)/24 + (J11*J12*J22*h^3)/6 + (J11*J12*h^2)/2 + (J21*J12^2*J22*h^4)/12 + (J21*J12^2*h^3)/6 + (J12*J22^3*h^4)/24 + (J12*J22^2*h^3)/6 + (J12*J22*h^2)/2 + J12*h)*linpoint[2]
  B1top = ((J11^3*h^4)/24 + (J11^2*h^3)/6 + (J12*J21*J11*h^4)/12 + (J11*h^2)/2 + (J12*J21*J22*h^4)/24 + (J12*J21*h^3)/6 + h)*B11
  B2top = ((J11^2*J12*h^4)/24 + (J11*J12*J22*h^4)/24 + (J11*J12*h^3)/6 + (J21*J12^2*h^4)/24 + (J12*J22^2*h^4)/24 + (J12*J22*h^3)/6 + (J12*h^2)/2)*B21
  B1 = B1top + B2top
  b1 = f0*h - D11*h - (D11*J11^2*h^3)/6 - (D11*J11^3*h^4)/24 + (J11^2*f0*h^3)/6 + (J11^3*f0*h^4)/24 - (D11*J11*h^2)/2 - (D21*J12*h^2)/2 + (J11*f0*h^2)/2 + (J12*g0*h^2)/2 + (J12*J21*f0*h^3)/6 + (J11*J12*g0*h^3)/6 + (J12*J22*g0*h^3)/6 - (D21*J11^2*J12*h^4)/24 - (D21*J12^2*J21*h^4)/24 - (D21*J12*J22^2*h^4)/24 + (J11^2*J12*g0*h^4)/24 + (J12^2*J21*g0*h^4)/24 + (J12*J22^2*g0*h^4)/24 - (D11*J12*J21*h^3)/6 - (D21*J11*J12*h^3)/6 - (D21*J12*J22*h^3)/6 - (D11*J11*J12*J21*h^4)/12 - (D11*J12*J21*J22*h^4)/24 - (D21*J11*J12*J22*h^4)/24 + (J11*J12*J21*f0*h^4)/12 + (J12*J21*J22*f0*h^4)/24 + (J11*J12*J22*g0*h^4)/24

  A21 = ((J11^3*J21*h^4)/24 + (J11^2*J21*J22*h^4)/24 + (J11^2*J21*h^3)/6 + (J12*J11*J21^2*h^4)/12 + (J11*J21*J22^2*h^4)/24 + (J11*J21*J22*h^3)/6 + (J11*J21*h^2)/2 + (J12*J21^2*J22*h^4)/12 + (J12*J21^2*h^3)/6 + (J21*J22^3*h^4)/24 + (J21*J22^2*h^3)/6 + (J21*J22*h^2)/2 + J21*h)*linpoint[1]
  A22 = ((J11^2*J12*J21*h^4)/24 + (J11*J12*J21*J22*h^4)/12 + (J11*J12*J21*h^3)/6 + (J12^2*J21^2*h^4)/24 + (J12*J21*J22^2*h^4)/8 + (J12*J21*J22*h^3)/3 + (J12*J21*h^2)/2 + (J22^4*h^4)/24 + (J22^3*h^3)/6 + (J22^2*h^2)/2 + J22*h)*linpoint[2]
  B1bot = ((J11^2*J21*h^4)/24 + (J11*J21*J22*h^4)/24 + (J11*J21*h^3)/6 + (J12*J21^2*h^4)/24 + (J21*J22^2*h^4)/24 + (J21*J22*h^3)/6 + (J21*h^2)/2)*B11
  B2bot = ((J22^3*h^4)/24 + (J22^2*h^3)/6 + (J12*J21*J22*h^4)/12 + (J22*h^2)/2 + (J11*J12*J21*h^4)/24 + (J12*J21*h^3)/6 + h)*B21
  B2 = B1bot + B2bot
  b2 = g0*h - D21*h - (D21*J22^2*h^3)/6 - (D21*J22^3*h^4)/24 + (J22^2*g0*h^3)/6 + (J22^3*g0*h^4)/24 - (D11*J21*h^2)/2 - (D21*J22*h^2)/2 + (J21*f0*h^2)/2 + (J22*g0*h^2)/2 + (J11*J21*f0*h^3)/6 + (J21*J22*f0*h^3)/6 + (J12*J21*g0*h^3)/6 - (D11*J11^2*J21*h^4)/24 - (D11*J12*J21^2*h^4)/24 - (D11*J21*J22^2*h^4)/24 + (J11^2*J21*f0*h^4)/24 + (J12*J21^2*f0*h^4)/24 + (J21*J22^2*f0*h^4)/24 - (D11*J11*J21*h^3)/6 - (D11*J21*J22*h^3)/6 - (D21*J12*J21*h^3)/6 - (D11*J11*J21*J22*h^4)/24 - (D21*J11*J12*J21*h^4)/24 - (D21*J12*J21*J22*h^4)/12 + (J11*J21*J22*f0*h^4)/24 + (J11*J12*J21*g0*h^4)/24 + (J12*J21*J22*g0*h^4)/12

  return [A11 A12;A21 A22], [B1;B2], [b1;b2]
end

end # module
