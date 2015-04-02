# This module contains the functions required to model
# the movements of the burglar in the house. Based on the problem in:
# Bayesian Reasoning and Machine Learning
# by David Barber
# Website: http://www0.cs.ucl.ac.uk/staff/d.barber/brml/
# Chapter 23, Example 23.3: Localisation Problem

module Burglar

using HMM

type house
  floor :: Array{Int64, 2}
  creaks :: Array{Int64, 2}
  bumps :: Array{Int64, 2}
  n :: Int64
  function house(n::Int64)
    f = zeros(Int64,n,n)
    c = rand(0:1, n,n)
    b = rand(0:1, n,n)
    initial_loc = rand(1:n)
    f[initial_loc] = 1
    new(f, c, b, n)
  end
end

function createhmm(model::house)
  tp = zeros(model.n^2, model.n^2)
  ep = zeros(4, model.n^2)

  # Transmission probability tables
  for col=1:model.n^2
    for row=1:model.n^2
      (row in getLegalMoves(model.n, col)) && (tp[row, col] += 1.0)
    end
  end
  tp = tp./sum(tp,1)

  # Emission probability tables
  cl = 0.9 # creak if light
  ncl = 1.0-cl
  cd = 0.01 # creak if dark
  ncd = 1.0-cd
  bl = 0.9 # bump if light
  nbl = 1.0-bl
  bd = 0.01 # bump if dark
  nbd = 1.0-bd
  # Emission 1: creak and bump
  # Emission 2: creak but no bump
  # Emission 3: no creak but a bump
  # Emission 4: neither a creak nor a bump
  for col=1:model.n^2
      if model.creaks[col] == 1 && model.bumps[col] == 1 # light light
        ep[:, col] = [cl*bl, cl*nbl, ncl*bl, ncl*nbl]
      elseif model.creaks[col] == 1 && model.bumps[col] == 0 # light dark
        ep[:, col] = [cl*bd, cl*nbd, ncl*bd, ncl*nbd]
      elseif model.creaks[col] == 0 && model.bumps[col] == 1 # dark light
        ep[:, col] = [cd*bl, cd*nbl, ncd*bl, ncd*nbl]
      else # dark dark
        ep[:, col] = [cd*bd, cd*nbd, ncd*bd, ncd*nbd]
      end
  end
  return HMM.hmm(tp, ep)
end

function move!(model::house)
  # Move the burglar
  moves = getLegalMoves(model)
  n = length(moves)
  newloc = moves[rand(1:n)]
  oldloc = getLocation(model)
  model.floor[oldloc] = 0
  model.floor[newloc] = 1

  return nothing
end

function getLegalMoves(model::house)
  # Returns an array of legal moves
  moves = Int64[]
  loc = getLocation(model)

  up = loc-1
  down = loc+1
  left = loc - model.n
  right = loc + model.n

  up in 1:model.n^2 && !(loc in 1:model.n:model.n^2) && push!(moves, up)
  down in 1:model.n^2 && !(loc in model.n:model.n:model.n^2) && push!(moves, down)
  left in 1:model.n^2 && push!(moves, left)
  right in 1:model.n^2 && push!(moves, right)

  return moves
end

function getLegalMoves(n::Int64, k::Int64)
  # Returns an array of legal moves
  moves = Int64[]
  loc = k

  up = loc-1
  down = loc+1
  left = loc - n
  right = loc + n

  up in 1:n^2 && !(loc in 1:n:n^2) && push!(moves, up)
  down in 1:n^2 && !(loc in n:n:n^2) && push!(moves, down)
  left in 1:n^2 && push!(moves, left)
  right in 1:n^2 && push!(moves, right)

  return moves
end

function getLocation(model::house)
  # Returns the location of the burglar
  return findfirst(model.floor, 1)
end

end #module
