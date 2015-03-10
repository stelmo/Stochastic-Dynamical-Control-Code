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
  nS, = size(model.A)

  for p=1:N
    for s=1:nS
      if particles.s[p] == s
        particles.w[p] = particles.w[p]*pdf(model.ydists[s], y - model.G[s](particles.x[:, p])) # weight of each particle
      end
    end
  end

  particles.w = particles.w ./ sum(particles.w)

  if numberEffectiveParticles(particles) < N/2
    resample!(particles)
  end
end

function filter!(particles::Particles, u, y, model::Model)

  nX, N = size(particles.x)
  nS, = size(model.A)

  # This can be made more compact but at the cost of clarity
  # first draw switch sample
  for p=1:N
    for s=1:nS
      if particles.s[p] == s
        particles.s[p] = rand(Categorical(model.A[:,s]))
      end
    end
  end

  # Now draw (predict) state sample
  for p=1:N
    for s=1:nS
      if particles.s[p] == s
        noise = rand(model.xdists[s])
        particles.x[:, p] = model.F[s](particles.x[:, p], u, noise) # predict
        particles.w[p] = particles.w[p]*pdf(model.ydists[s], y - model.G[s](particles.x[:, p])) # weight of each particle
      end
    end
  end

  particles.w = particles.w ./ sum(particles.w)
  (true in isnan(particles.w)) && warn("Particles have become degenerate! All is lost!")


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
    D = maximum(particles.x[k,:]) - minimum(particles.x[k,:])
    if D == 0.0
      warn("Particle distance very small! Roughening could cause problems...")
    end
    sig[k] = K*D*N^(-1./xN)
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
