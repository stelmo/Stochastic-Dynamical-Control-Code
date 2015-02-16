# Hidden Markov Model
# h - h -> ...
# |   |
# o   o
# Tests for the HMM code.
# All tests are compared to code supplied by
# Bayesian Reasoning and Machine Learning
# by David Barber
# Website: http://www0.cs.ucl.ac.uk/staff/d.barber/brml/

importall HMM_functions

A = [0.5 0.0 0.0;0.3 0.6 0.0;0.2 0.4 1.0]
B = [0.7 0.4 0.8;0.3 0.6 0.2]

mod1 = HMM(A, B)
initial = [0.9, 0.1, 0.0]
evidence = [1,2,1,1,1,1,1,1,1,1]

a1 = forward(mod1, initial, evidence) # works!

bw = zeros(length(initial), length(evidence))
for k=1:length(evidence)
  bw[:, k] = baum_welch(mod1, initial, evidence, k) # works!
end

vtb = viterbi(mod1, initial, evidence) # works!

pstates, pevidence = prediction(mod1, initial, evidence) # probably works!
