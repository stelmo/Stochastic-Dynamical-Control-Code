# Hidden Markov Model
# h - h -> ...
# |   |
# o   o
# Problem taken from:
# Bayesian Reasoning and Machine Learning
# by David Barber
# Website: http://www0.cs.ucl.ac.uk/staff/d.barber/brml/
# Chapter 23, Example 23.3: Localisation Problem

using Distributions
using PyPlot
importall HMM_functions
importall Burglar_functions

T = 10
n = 5

house = House(n)
hmm = createHMM(house) # create the hidden markov model

movements = zeros(Int64,n,n,T)
locs = zeros(Int64, T)
observations = zeros(Int64, T)

# Initialise
# initial = normalise(rand(n*n)) # no idea where the burglar is
initial = zeros(n*n)
initial[1:2*n] = 1.0/(2.*n)

# Measurement and actual movements
for t=1:T
  locs[t] = getLocation(house)
  movements[:,:,t] = house.floor
  observations[t] = rand(Categorical(hmm.ep[:,locs[t]]))
  move!(house)
end

## Inference
filter = forward(hmm, initial, observations)
bw = zeros(length(initial), length(observations))
for k=1:length(observations)
  bw[:, k] = baum_welch(hmm, initial, observations, k) # works!
end
vtb = viterbi(hmm, initial, observations)
mlmove = zeros(Int64,n,n,T)
for k=1:length(observations)
  temp = zeros(Int64, n,n)
  temp[vtb[k]] = 1.0
  mlmove[:,:, k] = temp
end

# Plotting
figure(1)
for t=1:T
  subplot(4, T, t)
  imshow(movements[:,:, t], cmap="Greys", interpolation="nearest")
  axis("off")
  subplot(4, T, t+T)
  imshow(reshape(filter[:, t], n,n), cmap="Greys",interpolation="nearest")
  axis("off")
  subplot(4, T, t+2*T)
  imshow(reshape(bw[:, t], n,n), cmap="Greys",interpolation="nearest")
  axis("off")
  subplot(4, T, t+3*T)
  imshow(mlmove[:,:, t], cmap="Greys",interpolation="nearest")
  axis("off")
end
