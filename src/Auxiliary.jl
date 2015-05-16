module Auxiliary

using Distributions
using KernelDensity
using PyPlot

function KL(part_states, part_weights, m, S)
  # Discrete Kullback-Leibler divergence test wrt a multivariate normal model.

  dnorm = MvNormal(m, S)
  dcat = Categorical(part_weights)
  N = length(part_weights)
  new_states = zeros(2, N)
  for k=1:N
    i = rand(dcat)
    new_states[:, k] = part_states[:, i]
  end
  estden = kde(new_states')

  kldiv = 0.0
  for k=1:N
    #draw from samplesw
    kldiv += -log(pdf(dnorm, new_states[:, k])) + log(pdf(estden, new_states[1, k], new_states[2, k]))
  end

  return (1.0/N)*kldiv
end

function showEstimatedDensity(part_states, part_weights)
  # Discrete Kullback-Leibler divergence test wrt a multivariate normal model.

  dcat = Categorical(part_weights)
  N = length(part_weights)
  new_states = zeros(2, N)
  for k=1:N
    i = rand(dcat)
    new_states[:, k] = part_states[:, i]
  end
  estden = kde(new_states')
  figure() # new figure otherwise very cluttered
  contour(estden)
end

end # module
