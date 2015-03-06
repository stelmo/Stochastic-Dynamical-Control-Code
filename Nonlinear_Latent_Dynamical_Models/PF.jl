module PF
# Particle Filter module as explained in the paper "Novel approach to nonlinear/
# non-Gaussian Bayesian state estimation" by Gorden et al (1993).

using Distributions

type Model
  f::Function # state transition model
  g::Function # state observation model
end

type Particle
  # Implements a particle
  x::Array{Float64, 1} # state
end

function init_PF(dist, nP::Int64, nY::Int64)
  # Initialise the particle filter.
  # dist => a distribution from the package Distributions. Implicitly specifies
  # the number of states per particle.
  # nY => number of observations per particle.
  # Return an array of nP particles.

  particles = Array(Particle, nP)

  for p=1:nP
    xs = rand(dist) # draw from the proposed prior
    ys = zeros(nY) # initialise to zero the observations
    w = 1.0/nP # uniform weight
    particles[p] = Particle(xs, ys, w)
  end

  return particles
end

function init_filter(particles::Array{Particle, 1}, u, y, plantdist, measuredist, model::Model)
  # Performs only the update step.
  N::Int64 = length(particles)
  w = zeros(N)

  for p=1:N
    w[p] = logpdf(measuredist, observation - model.g(particles[p].x)) # weight of each particle
  end

  w = w./cumsum(w) # normalise weights
  resample = rand(Categorical(w), N) # draw N samples from weighted Categorical
  copyparticles = copy(particles)
  for p=1:N # resample
    particles[p].x = copyparticles[resample[p]].x
  end

end

function filter(particles::Array{Particle, 1}, u, y, plantdist, measuredist, model::Model)
  # Performs the state prediction step.
  # plantnoise => distribution from whence the noise cometh
  N::Int64 = length(particles)
  w = zeros(N)

  for p=1:N
    noise = rand(plantdist)
    particles[p].x = model.f(particles[p].x, u, noise)
    w[p] = logpdf(measuredist, observation - model.g(particles[p].x)) # weight of each particle
  end

  w = w./cumsum(w) # normalise weights
  resample = rand(Categorical(w), N) # draw N samples from weighted Categorical
  copyparticles = copy(particles)
  for p=1:N # resample
    particles[p].x = copyparticles[resample[p]].x
  end
end

function roughen()

end

end
