# Hidden Markov Model
# h - h -> ...
# |   |
# o   o
# Tests for the HMM code.
# All tests are compared to code supplied by:
# Bayesian Reasoning and Machine Learning
# by David Barber
# Website: http://www0.cs.ucl.ac.uk/staff/d.barber/brml/
# The example is taken from Exercise 23.3
# The functions were all compared to their Matlab equivalent.
# I used: HMMforward, HMMsmooth and HMMviterbi

using Base.Test
importall HMM_functions

# Discrete model
A = [0.5 0.0 0.0;0.3 0.6 0.0;0.2 0.4 1.0] #transition probabilities
B = [0.7 0.4 0.8;0.3 0.6 0.2] # emission probabilities

mod1 = HMM(A, B) # create the HMM object
initial = [0.9, 0.1, 0.0] # initial state distribution
evidence = [1,1,2,1,2,1,2] # evidence/observations

filter_me = forward(mod1, initial, evidence)
filter_barber = readcsv("filter.csv") # read in the ideal answers

bw_me = zeros(length(initial), length(evidence))
for k=1:length(evidence)
  bw_me[:, k] = baum_welch(mod1, initial, evidence, k) # works!
end
bw_barber = readcsv("smooth.csv") #read in the ideal answers

vtb_me = viterbi(mod1, initial, evidence) # works!
vtb_barber = [1, 3, 3, 3, 3, 3, 3] # Barber's answer

pstates, pevidence = prediction(mod1, initial, evidence) # No test for this - not implemented by barber

# Run the tests
# Viterbi Inference
vtb_handler(r::Test.Success) = println("Successful Viterbi test!")
vtb_handler(r::Test.Failure) = error("Failure with the Viterbi test: $(r.expr)")
vtb_handler(r::Test.Error) = rethrow(r)
Test.with_handler(vtb_handler) do
  @test vtb_me == vtb_barber
end
# Filter Inference
filter_handler(r::Test.Success) = println("Successful filter test!")
filter_handler(r::Test.Failure) = error("Failure with the filter test: $(r.expr)")
filter_handler(r::Test.Error) = rethrow(r)
Test.with_handler(filter_handler) do
  @test maximum(abs(filter_me - filter_barber)) < 1e-4
end
# Smoother/Baum Welch Inference
smooth_handler(r::Test.Success) = println("Successful Baum-Welch test!")
smooth_handler(r::Test.Failure) = error("Failure with the Baum-Welch test: $(r.expr)")
smooth_handler(r::Test.Error) = rethrow(r)
Test.with_handler(smooth_handler) do
  @test maximum(abs(bw_me - bw_barber)) < 1e-4
end
