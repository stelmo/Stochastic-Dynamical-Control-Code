# Du Reactor
# Taken from the paper "Modeling and control of a continuous stirred tank reactor
# based on mixed logical dynamical model" by J. Du, C. Song and P. Li (2007). See
# the folder: literature/control on dropbox
module LinearReactor_functions

type LinearReactor
  # Supplies the parameters of the reactor. They may change i.e. drift...
  V :: Float64
  R :: Float64
  CA0 :: Float64
  TA0 :: Float64
  dH :: Float64
  k0 :: Float64
  E :: Float64
  Cp :: Float64
  rho :: Float64
  F :: Float64
end

function run_reactor(xprev::Array{Float64, 1}, u::Float64, h::Float64, model::DuReactor)
  # Use Runga-Kutta method to solve for the next time step using the full
  # nonlinear model.
  k1 :: Array{Float64, 1} = reactor_ode(xprev, u, model)
  k2 :: Array{Float64, 1} = reactor_ode(xprev + 0.5*h.*k1, u, model)
  k3 :: Array{Float64, 1} = reactor_ode(xprev + 0.5*h.*k2, u, model)
  k4 :: Array{Float64, 1} = reactor_ode(xprev + h.*k3, u, model)
  xnow :: Array{Float64, 1} = xprev + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
  return xnow
end

function reactor_ode(xprev::Array{Float64, 1}, u::Float64, model::DuReactor)
  # Evaluate the ODE functions describing the reactor.
  xnow :: Array{Float64, 1} = zeros(2)
  xnow[1] = -model.phi*xprev[1]*kappa(xprev[2], model) + model.q*(model.x1f-xprev[1])
  xnow[2] = model.beta*model.phi*xprev[1]*kappa(xprev[2], model) - (model.q+model.delta)*xprev[2]+model.delta*u + model.q*model.x2f
  return xnow
end

function kappa(x::Float64, model::DuReactor)
  # Evaluates the kappa function in the ode describing the reactor.
  kappa_eval :: Float64 = exp(x/(1.0+x/model.lambda))
  return kappa_eval
end

end # module
