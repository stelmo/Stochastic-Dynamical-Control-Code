# This script will run the linear Du reactor.

using PyPlot
import DuReactor_functions
reload("DuReactor_functions.jl")

# Specify the system parameters
du = begin
  phi = 0.072
  q = 1.0
  beta = 8.0
  delta = 0.3
  lambda = 20.0
  x1f = 1.0
  x2f = 0.0
  DuReactor_functions.DuReactor(phi, q, beta, delta, lambda, x1f, x2f)
end
