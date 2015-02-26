# Du Reactor
# Taken from the paper "Modeling and control of a continuous stirred tank reactor
# based on mixed logical dynamical model" by J. Du, C. Song and P. Li (2007). See
# the folder: literature/control on dropbox
module DuReactor

immutable DuReactor
  # Supplies the parameters of the reactor.
  phi :: Float64
  q :: Float64
  beta :: Float64
  delta :: Float64
  lambda :: Float64
  x1f :: Float64
  x2f :: Float64
end

function run_reactor(model::DuReactor)
  # Use Runga-Kutta method to solve for the next time step using the full
  # nonlinear model.
end

end # module
