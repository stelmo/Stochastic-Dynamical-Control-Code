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
  w::Array{Float64, 1} # collection of particle weights
end

function init_PF(dist, nP::Int64, xN::Int64)
  # Initialise the particle filter.
  # dist => a distribution from the package Distributions.
  # nX => number of states per particle.
  # Return an array of nP particles.

  particles = Particles(zeros(xN, nP), zeros(nP))
  for p=1:nP
    drawx = rand(dist) # draw from the proposed prior
    particles.x[:, p] = drawx
    particles.w[p] = 1./nP # uniform initial weight
  end

  return particles
end

function init_filter!(particles::Particles, u, y, measuredist, model::Model)
  # Performs only the update step.
  nX, N = size(particles.x)

  for p=1:N
    particles.w[p] = particles.w[p]*pdf(measuredist, y - model.g(particles.x[:, p])) # weight of each particle
  end

  (abs(maximum(particles.w)) < 1e-8) && warn("The particles all have very small weight...")
  particles.w = particles.w ./ sum(particles.w)
  (true in isnan(particles.w)) && error("Particles have become degenerate!")

  if numberEffectiveParticles(particles) < N/2
    resample!(particles)
  end
end

function resample!(particles::Particles)
  N = length(particles.w)
  resample = rand(Categorical(particles.w), N) # draw N samples from weighted Categorical
  copyparticles = copy(particles.x)
  for p=1:N # resample
    particles.x[:,p] = copyparticles[:, resample[p]]
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

function filter!(particles::Particles, u, y, plantdist, measuredist, model::Model)
  # Performs the state prediction step.
  # plantnoise => distribution from whence the noise cometh
  # measuredist => distribution from whence the plant uncertainty cometh

  nX, N = size(particles.x)

  for p=1:N
    noise = rand(plantdist)
    particles.x[:, p] = model.f(particles.x[:, p], u, noise) # predict
    particles.w[p] = particles.w[p]*pdf(measuredist, y - model.g(particles.x[:, p])) # weight of each particle
  end

  (abs(maximum(particles.w)) < 1e-8) && warn("The particles all have very small weight...")
  particles.w = particles.w ./ sum(particles.w)
  (true in isnan(particles.w)) && error("Particles have become degenerate!")

  if numberEffectiveParticles(particles) < N/2
    resample!(particles)
  end
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

function predict!(parts, u, plantdist, model)
  # Project the particles one step forward.
  # NOTE: this overwrite parts therefore use a dummy variable!

  nX, nP = size(parts)
  for p=1:nP
      noise = rand(plantdist)
      parts[:, p] = model.f(parts[:, p], u, noise) # predict
  end

end

end # Module
