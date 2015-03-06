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
  y::Array{Float64, 1} # observation
  w::Float64 # weight
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

function predict!(particles::Array{Particle, 1}, model::Model, plantnoise)
  # Performs the state prediction step.
  # dist => distribution from whence the noise cometh
  N::Int64 = length(particles)
  for p=1:N
    noise = rand(plantnoise)
    particles[p].x = model.f(particles[p].x, noise)
  end
end

function update!(particles::Array{Particle, 1}, model::Model, measurementnoise, observation)
  # Performs the Bayesian update step given observations.
  
end

function roughen()

end



end
