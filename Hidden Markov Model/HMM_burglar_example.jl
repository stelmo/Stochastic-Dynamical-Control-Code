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
using HMM_functions
using Burglar_functions

T = 10
n = 5

house = House(n)
hmm = createHMM(house) # create the hidden markov model

movements = zeros(Int64,n,n,T)
locs = zeros(Int64, T)
observations = zeros(Int64, T)

# Initial distribution
# initial = normalise(rand(n*n)) # no idea where the burglar is - not good
initial = zeros(n*n) # know the burglar is near the entrance of the house
# initial[1:2*n] = 1.0/(2.*n) # two column initial guess
initial[1:n] = 1.0/(n) # one column initial guess

# Measurement and actual movements
for t=1:T
  locs[t] = getLocation(house)
  movements[:,:,t] = house.floor
  observations[t] = rand(Categorical(hmm.ep[:,locs[t]]))
  move!(house)
end

## Inference
# Forward Filter
filter = forward(hmm, initial, observations)

# Smoothing
fbs = zeros(length(initial), length(observations))
for k=1:length(observations)
  fbs[:, k] = smooth(hmm, initial, observations, k)
end

# Viterbi
vtb = viterbi(hmm, initial, observations)
mlmove = zeros(Int64,n,n,T) # construct floor matrices showing the viterbi path
for k=1:length(observations)
  temp = zeros(Int64, n,n)
  temp[vtb[k]] = 1.0
  mlmove[:,:, k] = temp
end

# Prediction
predmove = zeros(n, n, T)
predmove[:,:,1] = reshape(initial, n, n) # first time step is just the prior
for k=2:T
  pstate, = prediction(hmm, initial, observations[1:k-1])
  predmove[:,:, k] = round(reshape(pstate, n, n), 2) # round to make predictions stand out more
end

fs = 18 #font size
figure(1) # Inference - no prediction
for t=1:T
  subplot(4, T, t)
  imshow(movements[:,:, t], cmap="Greys", interpolation="nearest")
  title("t=$(t)",fontsize=fs)
  t==1 && ylabel("True Location", fontsize=fs)
  tick_params(axis="both", which="both", bottom="off", top="off", left="off", right="off", labelbottom="off", labelleft="off")

  subplot(4, T, t+T)
  imshow(reshape(filter[:, t], n,n), cmap="Greys",interpolation="nearest")
  t==1 && ylabel("Filtering", fontsize=fs)
  tick_params(axis="both", which="both", bottom="off", top="off", left="off", right="off", labelbottom="off", labelleft="off")

  subplot(4, T, t+2*T)
  imshow(reshape(fbs[:, t], n,n), cmap="Greys",interpolation="nearest")
  t==1 && ylabel("Smoothing", fontsize=fs)
  tick_params(axis="both", which="both", bottom="off", top="off", left="off", right="off", labelbottom="off", labelleft="off")

  subplot(4, T, t+3*T)
  imshow(mlmove[:,:, t], cmap="Greys",interpolation="nearest")
  t==1 && ylabel("Viterbi Inference", fontsize=fs)
  tick_params(axis="both", which="both", bottom="off", top="off", left="off", right="off", labelbottom="off", labelleft="off")
end

figure(2) # Inference - prediction
for t=1:T
  subplot(2, T, t)
  imshow(movements[:,:, t], cmap="Greys", interpolation="nearest")
  title("t=$(t)", fontsize=fs)
  t==1 && ylabel("True Location", fontsize=fs)
  tick_params(axis="both", which="both", bottom="off", top="off", left="off", right="off", labelbottom="off", labelleft="off")

  subplot(2, T, t+T)
  imshow(predmove[:,:, t], cmap="Greys", interpolation="nearest")
  t==1 && ylabel("Predicted Location", fontsize=fs)
  tick_params(axis="both", which="both", bottom="off", top="off", left="off", right="off", labelbottom="off", labelleft="off")
end
