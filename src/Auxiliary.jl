module Auxiliary

using Distributions
using KernelDensity
using PyPlot

warn("Auxiliary is hardcoded for the CSTR!")
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

  sweights = 1.0 - sum(part_weights)
  if sweights < 0.0 # add some robustness here (very crude)
    part_weights[1] = part_weights[1] + sweights
  elseif sweights > 0.0
    part_weights[1] = part_weights[1] - sweights
  end


  dcat = Categorical(part_weights)
  N = length(part_weights)
  new_states = zeros(2, N)
  for k=1:N
    i = rand(dcat)
    new_states[:, k] = part_states[:, i]
  end
  estden = kde(new_states')
  rc("font", family="serif", size=24)
  figure() # new figure otherwise very cluttered
  contour(estden)
end

end # module
