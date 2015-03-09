module SPF
# switching particle filter

type Particles
  x :: Array{Float64, 2} # states
  s :: Array{Int64, 1} # switches
  w :: Array{Float64, 1} # weights
end

type Model
  F :: Array{Function, 1} # transition
  G :: Array{Function, 1} # emission
  ydist # no idea what this type is
end

type Switch
  A :: Array{Float64, 2} # HMM model, columns sum to 1
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





end #module
