module Auxiliary

using Distributions
using KernelDensity
using PyPlot

warn("Auxiliary is hardcoded for the CSTR!")
function KL(part_states, part_weights, m, S, temp_states)
  # Discrete Kullback-Leibler divergence test wrt a multivariate normal model.

  sweights = 1.0 - sum(part_weights)
  N = length(part_weights)
  if sweights < 0.0 # add some robustness here (very crude)
    part_weights = part_weights .+ (1.0/N)*sweights
    w = string("Particle weights adjusted by ", sweights, " in Auxiliary!")
    warn(w)
  elseif sweights > 0.0
    part_weights = part_weights .+ (1.0/N)*sweights
    w = string("Particle weights adjusted by ", sweights, " in Auxiliary!")
    warn(w)
  end

  dnorm = MvNormal(m, S)
  dcat = Categorical(part_weights)

  for k=1:N
    i = rand(dcat)
    temp_states[:, k] = part_states[:, i]
  end
  estden = kde(temp_states')

  kldiv = 0.0
  for k=1:N
    #draw from samplesw
    kldiv += -log(pdf(dnorm, temp_states[:, k])) + log(pdf(estden, temp_states[1, k], temp_states[2, k]))
  end

  return (1.0/N)*kldiv
end

function KLbase(m, S, temp_states, N)


    dnorm = MvNormal(m, S)

    for k=1:N
      temp_states[:, k] = rand(dnorm)
    end
    estden = kde(temp_states')

    kldiv = 0.0
    for k=1:N
      #draw from samplesw
      kldiv += -log(pdf(dnorm, temp_states[:, k])) + log(pdf(estden, temp_states[1, k], temp_states[2, k]))
    end

    return (1.0/N)*kldiv
end

function KLuniform(m, S, temp_states, N)

    s11 = S[1]
    s22 = S[4]

    dnorm = MvNormal(m, S)

    m1 = [m[1]-sqrt(s11), m[1]+sqrt(s11)]
    m2 = [m[2]-sqrt(s22), m[2]+sqrt(s22)]

    dnorm1 = Uniform(minimum(m1)*2, maximum(m1)*2)
    dnorm2 = Uniform(minimum(m2)*2, maximum(m2)*2)

    for k=1:N
      temp_states[1, k] = rand(dnorm1)
      temp_states[2, k] = rand(dnorm2)
    end
    estden = kde(temp_states')

    kldiv = 0.0
    for k=1:N
      #draw from samplesw
      kldiv += -log(pdf(dnorm, temp_states[:, k])) + log(pdf(estden, temp_states[1, k], temp_states[2, k]))
    end

    return (1.0/N)*kldiv
end

function showEstimatedDensity(part_states, part_weights, temp_states)
  # Discrete Kullback-Leibler divergence test wrt a multivariate normal model.

  sweights = 1.0 - sum(part_weights)
  N = length(part_weights)
  if sweights < 0.0 # add some robustness here (very crude)
    part_weights = part_weights .+ (1.0/N)*sweights
    w = string("Particle weights adjusted by ", sweights, " in Auxiliary!")
    warn(w)
  elseif sweights > 0.0
    part_weights = part_weights .+ (1.0/N)*sweights
    w = string("Particle weights adjusted by ", sweights, " in Auxiliary!")
    warn(w)
  end

  dcat = Categorical(part_weights)

  for k=1:N
    i = rand(dcat)
    temp_states[:, k] = part_states[:, i]
  end
  estden = kde(temp_states')
  rc("text", usetex=true)
  rc("font", family="serif", serif="Computer Modern", size=24)

  figure() # new figure otherwise very cluttered
  contour(estden)
  xlabel(L"C_A [kmol.m^{-3}]")
  ylabel(L"T_R [K]")
end

function removeOutliers(xs, multiple=2)
  # remove columns if the column has an element which is more than twice as the mean

  rows, cols = size(xs)
  m1 = mean(xs, 2)
  remindex = Int64[]
  for k=1:cols
    if xs[1, k] > m1[1]*multiple || xs[2, k] > m1[2]*multiple
      push!(remindex, k)
    end
  end

  fxs = zeros(rows, cols-length(remindex))
  counter = 1
  for k=1:cols
    if !(k in remindex)
    fxs[:, counter] = xs[:, k]
    counter += 1
    end
  end

  return fxs
end

end # module
