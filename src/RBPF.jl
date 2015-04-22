# Rao Blackwellised Particle Filter
module RBPF
# WARNING: this is made specifically for the system I am investigating
warn("RBPF is hardcoded for this system!")
using Distributions
using Reactor
using LLDS
using SPF

type Particles
  mus :: Array{Float64, 2} # mean
  sigmas :: Array{Float64, 3} # covariance
  ss :: Array{Int64, 1} # switches
  ws :: Array{Float64, 1} # weights
end

type Model
  A
  B
  b
  C
  Q
  R
end

function setup_RBPF(linsystems, C, Q, R)
  # Setup each switch
  N = length(linsystems)
  models = Array(Model, N)
  for k=1:N
    models[k] = Model(linsystems[k].A, linsystems[k].B, linsystems[k].b, C, Q, R)
  end
  A = SPF.calcA(linsystems)
  return models, A
end

function init_RBPF(sdist, mu_init, sigma_init, xN, nP)
  # Initialise the particle filter.

  particles = Particles(zeros(xN, nP), zeros(xN, xN, nP), zeros(Int64, nP), zeros(nP))
  for p=1:nP
    sdraw = rand(sdist)
    particles.mus[:, p] = mu_init # normal mu
    particles.sigmas[:,:, p] = sigma_init
    particles.ss[p] = sdraw
    particles.ws[p] = 1./nP # uniform initial weight
  end

  return particles
end

function init_filter!(particles::Particles, u, y, models::Array{Model, 1})

  nX, N = size(particles.mus)
  nS = length(models)

  for p=1:N
    for s=1:nS
      if particles.ss[p] == s
        particles.mus[:,p] = particles.mus[:, p] - models[particles.ss[p]].b # adjust mu for specific switch
        mu = models[s].C * (models[s].A*particles.mus[:, p] + models[s].B*u)
        sigma = models[s].C * (models[s].A*particles.sigmas[:,:,p]*models[s].A' + models[s].Q) * models[s].C' + models[s].R
        d = MvNormal(mu, sigma)
        if length(y) == 1
          particles.ws[p] = particles.ws[p]*pdf(d, [y-models[s].b[2]]) # HARDCODED for this system!!! ****************************NB
        else
          particles.ws[p] = particles.ws[p]*pdf(d, y-models[s].b) # weight of each particle
        end
        # println("Switch: ", s, " Predicts: ", round(mu + models[s].b, 4), "Observed: ", round(y,4), " Weight: ", round(particles.ws[p], 5))
        particles.mus[:,p] = particles.mus[:, p] + models[particles.ss[p]].b # fix mu for specific switch

      end
    end
  end


  particles.ws = particles.ws ./ sum(particles.ws)

  if numberEffectiveParticles(particles) < N/2
    resample!(particles, models)
  end
end

function filter!(particles::Particles, u, y, models::Array{Model, 1}, A)

  nX, N = size(particles.mus)
  nS = length(models)

  # This can be made more compact but at the cost of clarity
  # first draw switch sample
  for p=1:N
    for s=1:nS
      if particles.ss[p] == s
        particles.ss[p] = rand(Categorical(A[:,s]))
      end
    end
  end

  # apply KF and weight
  for p=1:N
    for s=1:nS
      if particles.ss[p] == s
        mu = models[s].C * (models[s].A*(particles.mus[:, p] - models[s].b) + models[s].B*u)
        sigma = models[s].C * (models[s].A*particles.sigmas[:,:,p]*models[s].A' + models[s].Q) * models[s].C' + models[s].R
        d = MvNormal(mu, sigma)
        if length(y) == 1
          particles.ws[p] = particles.ws[p]*pdf(d, [y-models[s].b[2]]) # HARDCODED for this system!!! ****************************NB
        else
          particles.ws[p] = particles.ws[p]*pdf(d, y-models[s].b) # weight of each particle
        end
        (isnan(particles.ws[p])) && (warn("Particle weight issue..."); particles.ws[p] = 0.0) # in case

        # println("Switch: ", s, " Predicts: ", round(mu + models[s].b, 4), "Observed: ", round(y,4), " Weight: ", round(particles.ws[p], 5))


        pmean = models[s].A*(particles.mus[:,p] - models[s].b) + models[s].B*u
        pvar =  models[s].Q + models[s].A*particles.sigmas[:,:, p]*transpose(models[s].A)
        kalmanGain = pvar*transpose(models[s].C)*inv(models[s].C*pvar*transpose(models[s].C) + models[s].R)
        ypred = models[s].C*pmean #predicted measurement
        if length(y) == 1
          updatedMean = pmean + kalmanGain*(y - models[s].b[2] - ypred) # adjust for state space
        else
          updatedMean = pmean + kalmanGain*(y - models[s].b - ypred) # adjust for state space
        end
        rows, cols = size(pvar)
        updatedVar = (eye(rows) - kalmanGain*models[s].C)*pvar

        particles.sigmas[:,:, p] = updatedVar
        particles.mus[:, p] = updatedMean + models[s].b # fix
      end
    end
  end

  # particles.ws = particles.ws .+ abs(minimum(particles.ws)) #no negative number issue
  (maximum(particles.ws) < (1.0/(N^2))) && warn("The particles all have very small weight...")
  particles.ws = particles.ws ./ sum(particles.ws)
  (true in isnan(particles.ws)) && error("Particles have become degenerate!")


  if numberEffectiveParticles(particles) < N/2
    resample!(particles,models)
  end
end

function resample!(particles::Particles, models::Array{Model, 1})
  N = length(particles.ws)
  resample = rand(Categorical(particles.ws), N) # draw N samples from weighted Categorical
  copyparticles_mus = copy(particles.mus)
  copyparticles_sigmas = copy(particles.sigmas)
  copyparticles_ss = copy(particles.ss)
  for p=1:N # resample
    particles.mus[:,p] = copyparticles_mus[:, resample[p]]
    particles.sigmas[:,:, p] = copyparticles_sigmas[:,:, resample[p]]
    particles.ss[p] = copyparticles_ss[resample[p]]
    particles.ws[p] = 1./N
  end
  roughen!(particles)
end

function numberEffectiveParticles(particles::Particles)
  # Return the effective number of particles.
  N = length(particles.ws)
  numeff = 0.0
  for p=1:N
    numeff += particles.ws[p]^2
  end
  return 1./numeff
end

function roughen!(particles::Particles)
  # Roughening the samples to promote diversity
  xN, N = size(particles.mus)
  sig = zeros(xN)

  K= 0.2 # parameter...
  flag = true
  for k=1:xN
    D = maximum(particles.mus[k,:]) - minimum(particles.mus[k,:])
    if D == 0.0
      warn("Particle distance very small! Roughening could cause problems...")
      flag = false
      break
    end
    sig[k] = K*D*N^(-1./xN)
  end

  if flag
    sigma = diagm(sig.^2)
    jitter = MvNormal(sigma)
    for p=1:N
      particles.mus[:, p] = particles.mus[:, p] + rand(jitter)
    end
  end
end

function getAveStats(particles::Particles)
  nX, nP = size(particles.mus)
  ave = zeros(nX)
  avesigma = zeros(nX, nX)
  for p=1:nP
    ave = ave + particles.ws[p]*particles.mus[:, p]
    avesigma = avesigma + particles.ws[p].*particles.sigmas[:,:, p]
  end

  return ave, avesigma
end

function getMLStats(particles::Particles)
  nX, nP = size(particles.mus)
  mlmu = zeros(nX)
  mlsigma = zeros(nX, nX)
  prevmaxweight = 0.0
  for p=1:nP
    if particles.ws[p] > prevmaxweight
      mlmu = particles.mus[:, p]
      mlsigma = particles.sigmas[:,:, p]
      prevmaxweight = particles.ws[p]
    end
  end

  return mlmu, mlsigma
end

function getMaxTrack(particles, numSwitches)
  maxtrack = zeros(numSwitches)
  numParticles = length(particles.ws)
  totals = zeros(numSwitches)
  for p=1:numParticles
    totals[particles.ss[p]] += particles.ws[p]
  end

  maxtrack[indmax(totals)] = 1.0
  return maxtrack
end

function smoothedTrack(numSwitches, switchtrack, ind, N)
  # returns a smoothed version of maxtrack given the history of the switch movements
  sN = getHistory(ind, N) # number of time intervals backwards we look
  modswitchtrack = sum(switchtrack[:, ind-sN:ind], 2)
  modmaxtrack = zeros(numSwitches)
  modmaxtrack[indmax(modswitchtrack)] = 1.0
  return modmaxtrack
end

function getHistory(ind, N)
  history = 0
  flag = true

  for k=1:N
    if ind-k == 0
      history = k
      flag = false
      break
    end
  end

  if flag
    history = N
  end

  return history-1
end

function getInitialSwitches(initial_states, linsystems::Array{Reactor.LinearReactor,1})
  N = length(linsystems)
  initstates = zeros(N) #pre-allocate

  for i=1:N
    initstates[i] = norm(linsystems[i].op-initial_states)
  end

  a = sort(initstates, rev=true)
  posA = zeros(N)
  for i=1:N
    posA[i] = find((x)->x==initstates[i], a)[1]
  end

  posA = posA./sum(posA,1)

  return posA
end

end #module
