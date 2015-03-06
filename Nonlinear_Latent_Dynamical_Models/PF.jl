module PF
# Particle Filter module as explained in the paper "Novel approach to nonlinear/
# non-Gaussian Bayesian state estimation" by Gorden et al (1993).

using Distributions

type Model
  f::Function # state transition model
  g::Function # state observation model
end

type Particles
  # Implements a particle
  x::Array{Float64, 2} # collection of particles
end

function init_PF(dist, nP::Int64, xN::Int64)
  # Initialise the particle filter.
  # dist => a distribution from the package Distributions.
  # nX => number of states per particle.
  # Return an array of nP particles.

  particles = Particles(zeros(xN, nP))
  for p=1:nP
    drawx = rand(dist) # draw from the proposed prior
    particles.x[:, p] = drawx
  end

  return particles
end

function init_filter!(particles::Particles, u, y, plantdist, measuredist, model::Model)
  # Performs only the update step.
  nX, N = size(particles.x)
  w = zeros(N)

  for p=1:N
    w[p] = pdf(measuredist, y - model.g(particles.x[:, p])) # weight of each particle
  end

  w = w./sum(w) # normalise weights
  resample = rand(Categorical(w), N) # draw N samples from weighted Categorical
  copyparticles = copy(particles.x)
  for p=1:N # resample
    particles.x[:,p] = copyparticles[:, resample[p]]
  end
end

function filter!(particles::Particles, u, y, plantdist, measuredist, model::Model)
  # Performs the state prediction step.
  # plantnoise => distribution from whence the noise cometh
  nX, N = size(particles.x)
  w = zeros(N)

  for p=1:N
    noise = rand(plantdist)
    particles.x[:, p] = model.f(particles.x[:, p], u, noise) # predict
    w[p] = pdf(measuredist, y - model.g(particles.x[:, p])) # update weight of each particle
  end

  w = w./sum(w) # normalise weights
  resample = rand(Categorical(w), N) # draw N samples from weighted Categorical
  copyparticles = copy(particles.x)
  for p=1:N # resample
    particles.x[:,p] = copyparticles[:, resample[p]]
  end
end

function getStats(particles::Particles)

  fitted = fit(MvNormal, particles.x)

  return mean(fitted), cov(fitted)
end


end # Module
