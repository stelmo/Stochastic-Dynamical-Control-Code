# Particle Swarm Optimisation
module PSO

using PF
using Distributions

function fitness(particles, us, ysp, usp, Q, R, plantdist, model)
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

  for n=1:N # predicted error
    PF.predict!(parts, us[n], plantdist, model) # modifies parts
    for p=1:nP
      fit = fit + particles.w[p]*((parts[:, p]-ysp)'*Q*(parts[:, p]-ysp)) + (us[n]-usp)'*R*(us[n]-usp)
    end
  end

  return fit
end



end # PSO module
