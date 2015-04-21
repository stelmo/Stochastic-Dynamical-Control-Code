# Particle Swarm Optimisation
module PSO

using PF
using Distributions
using Reactor
using NLsolve

type Particle
  pos
  vel
  prevbest
  prevbestpos
end

type Solution
  pos
  posfit
end


function fitness(particles, us, ysp, usp, Q, R, plantdist, model, skip)
  # This should ideally be run a few times to get the MC expected value.
  # NOTE: we use the standard quadratic function in the expectation i.e.
  # E[(x-y)'Q(x-y)]

  parts = copy(particles.x)
  fit = 0.0
  N = length(us)
  nX, nP = size(parts)

  for p=1:nP # initial error
    fit = fit + particles.w[p]*((parts[:, p]-ysp)'*Q*(parts[:, p]-ysp))
  end

  mvcounter = 1
  for n=1:N*(1+skip)-1 # predicted error

    #update mvcounter
    if n%(1+skip) == 0
      mvcounter += 1
    end

    PF.predict!(parts, us[mvcounter], plantdist, model) # modifies parts
    for p=1:nP
      fit = fit + particles.w[p]*(parts[:, p]-ysp)'*Q*(parts[:, p]-ysp)
    end

    fit = fit + (us[mvcounter]-usp)'*R*(us[mvcounter]-usp)
  end

  return fit/N
end

function parFitness(particles, us, ysp, usp, Q, R, plantdist, reactormodel, h)
  # A parallel implementation of the fitness function - suprisingly not so fast...
  nX, nP = size(particles.x)
  parts = copy(particles.x)
  fit = 0.0
  N = length(us)

  for p=1:nP # initial error
    fit = fit + particles.w[p]*((parts[p]-ysp)'*Q*(parts[p]-ysp))
  end

  for n=1:N # predicted error
    PF.parPredict!(parts, us[n], plantdist, reactormodel, h) # modifies parts
    for p=1:nP
      fit = fit + particles.w[p]*((parts[p]-ysp)'*Q*(parts[p]-ysp))
    end
    fit = fit + (us[mvcounter]-usp)'*R*(us[mvcounter]-usp)

  end

  return fit/N
end

function offset(y_concentration, model)
  # Finds the steady state given the controller input. Assume only concentration may be controlled.
  f!(x, fvec) = Reactor.reactor_func!([y_concentration, x[1]], x[2], model, fvec)
  offset_results = nlsolve(f!, [0.5, 0.0])
  return offset_results
end

function createswarm(nP, nX, umin, umax)
  # Initialise the swarm with random positions and velocities within a range.
  swarm = Array(Particle, nP) # create swarm
  for p=1:nP
    pnow = randomr(nX, umin, umax)
    swarm[p] = Particle(pnow, randomr(nX, umin, umax), 0.0, pnow)
  end
  sol = Solution(randomr(nX, umin, umax), 0.0)
  return swarm, sol
end

function initswarm(nP, nX, umin, umax, particles, ysp, usp, Q, R, plantdist, model, skip)

  swarm = Array(Particle, nP) # create swarm
  sol = Solution(randomr(nX, umin, umax), Inf)

  for p=1:nP
    pnow = randomr(nX, umin, umax)
    swarm[p] = Particle(pnow, randomr(nX, umin, umax), 0.0, pnow)
    swarm[p].prevbest = fitness(particles, swarm[p].pos, ysp, usp, Q, R, plantdist, model, skip)
    if swarm[p].prevbest[1] < sol.posfit[1]
      sol.pos = swarm[p].prevbestpos
      sol.posfit = swarm[p].prevbest
    end
  end
  return swarm, sol
end

function updateswarm!(swarm, sol, particles, ysp, usp, Q, R, plantdist, model, skip)
# Using optimum parameters
  nP = length(swarm)
  nX = length(swarm[1].pos)
  for p=1:nP
    swarm[p].vel = 0.7298.*swarm[p].vel + randomr(nX, 0.0, 1.49618).*(swarm[p].prevbestpos - swarm[p].pos) +  randomr(nX, 0.0, 1.49618).*(sol.pos - swarm[p].pos)
    swarm[p].pos = swarm[p].pos + swarm[p].vel

    tempfitness = fitness(particles, swarm[p].pos, ysp, usp, Q, R, plantdist, model, skip)

    if tempfitness[1] < sol.posfit[1]
      sol.pos = swarm[p].pos
      sol.posfit = tempfitness
    end

    if tempfitness[1] < swarm[p].prevbest[1]
      swarm[p].prevbestpos = swarm[p].pos
      swarm[p].prevbest = tempfitness
    end

  end

end

function optimise!(swarm, sol, particles, ysp, usp, Q, R, plantdist,model, skip)
  # Run the PSO
  for k=1:100
    updateswarm!(swarm, sol, particles, ysp, usp, Q, R, plantdist, model, skip)
    # println(sol.posfit[1])
  end
  # println("***********")
  return sol.pos[1]
end

function randomr(nX, xmin, xmax)
  # Returns a uniformly distributed random number in the range (xmin, xmax]
  return rand(nX).*(xmax-xmin) .+ xmin
end


end # PSO module
