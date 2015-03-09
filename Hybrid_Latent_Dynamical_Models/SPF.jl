module SPF
# switching particle filter

using Distributions

type Particles
  x :: Array{Float64, 2} # states
  s :: Array{Int64, 1} # switches
  w :: Array{Float64, 1} # weights
end

type Model
  F :: Array{Function, 1} # transition
  G :: Array{Function, 1} # emission
  A :: Array{Float64, 2} # HMM model, columns sum to 1
  xdists
  ydists # no idea what this type is
end

function init_SPF(xdist, sdist, nP::Int64, xN::Int64)
  # Initialise the particle filter.
  # xdist => state prior (a distribution from the package Distributions)
  # sdist => switch prior distribution
  # nP => number of particles
  # nX => number of states per particle
  # Return an array of nP particles

  particles = Particles(zeros(xN, nP), zeros(Int64, nP), zeros(nP))
  for p=1:nP
    xdraw = rand(xdist)
    sdraw = rand(sdist)
    particles.x[:, p] = xdraw
    particles.s[p] = sdraw
    particles.w[p] = 1./nP # uniform initial weight
  end

  return particles
end

function init_filter!(particles::Particles, u, y, model::Model)

  nX, N = size(particles.x)

  for p=1:N
    if particles.s[p] == 1
      particles.w[p] = particles.w[p]*pdf(model.ydists[1], y - model.G[1](particles.x[:, p])) # weight of each particle
    elseif particles.s[p] == 2
      particles.w[p] = particles.w[p]*pdf(model.ydists[2], y - model.G[2](particles.x[:, p])) # weight of each particle
    else
      particles.w[p] = particles.w[p]*pdf(model.ydists[3], y - model.G[3](particles.x[:, p])) # weight of each particle
    end
  end

  particles.w = particles.w ./ sum(particles.w) # normalise weights

  if numberEffectiveParticles(particles) < N/2
    resample!(particles)
  end
end

function filter!(particles::Particles, u, y, model::Model)

  nX, N = size(particles.x)

  # This can be made more compact but at the cost of clarity
  # first draw switch sample
  for p=1:N
    if particles.s[p] == 1
      particles.s[p] = rand(Categorical(model.A[:,1]))
    elseif particles.s[p] == 2
      particles.s[p] = rand(Categorical(model.A[:,2]))
    else
      particles.s[p] = rand(Categorical(model.A[:,3]))
    end
  end

  # Now draw (predict) state sample
  for p=1:N
    if particles.s[p] == 1
      noise = rand(model.xdists[1])
      particles.x[:, p] = model.F[1](particles.x[:, p], u, noise) # predict
      particles.w[p] = particles.w[p]*pdf(model.ydists[1], y - model.G[1](particles.x[:, p])) # weight of each particle
    elseif particles.s[p] == 2
      noise = rand(model.xdists[2])
      particles.x[:, p] = model.F[2](particles.x[:, p], u, noise) # predict
      particles.w[p] = particles.w[p]*pdf(model.ydists[2], y - model.G[2](particles.x[:, p])) # weight of each particle
    else
      noise = rand(model.xdists[3])
      particles.x[:, p] = model.F[3](particles.x[:, p], u, noise) # predict
      particles.w[p] = particles.w[p]*pdf(model.ydists[3], y - model.G[3](particles.x[:, p])) # weight of each particle
    end
  end

  particles.w = particles.w ./ sum(particles.w) # normalise weights

  if numberEffectiveParticles(particles) < N/2
    resample!(particles)
  end
end


function resample!(particles::Particles)
  N = length(particles.w)
  resample = rand(Categorical(particles.w), N) # draw N samples from weighted Categorical
  copyparticles_x = copy(particles.x)
  copyparticles_s = copy(particles.s)
  for p=1:N # resample
    particles.x[:,p] = copyparticles_x[:, resample[p]]
    particles.s[p] = copyparticles_s[resample[p]]
    particles.w[p] = 1./N
  end
  roughen!(particles)
end

function numberEffectiveParticles(particles::Particles)
  # Return the effective number of particles.
  N = length(particles.w)
  numeff = 0.0
  for p=1:N
    numeff += particles.w[p]^2
  end
  return 1./numeff
end

function roughen!(particles::Particles)
  # Roughening the samples to promote diversity
  xN, N = size(particles.x)
  sig = zeros(xN)

  K= 0.2 # parameter...

  for k=1:xN
    sig[k] = K*(maximum(particles.x[k,:]) - minimum(particles.x[k,:]))*N^(-1./xN)
  end

  sigma = diagm(sig.^2)
  jitter = MvNormal(sigma)
  for p=1:N
    particles.x[:, p] = particles.x[:, p] + rand(jitter)
  end
end

function getStats(particles::Particles)
  # Return the Gaussian statistics of the particles.
  fitted = fit(MvNormal, particles.x, particles.w)
  return mean(fitted), cov(fitted)
end

end #module
