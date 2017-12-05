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

if VERSION >= v"0.5-"
    using BaseTestDeprecated
    const Test = BaseTestDeprecated
else
    using Base.Test
end

using HMM

# First get current working directory - to make Travis happy
dircontent = readdir()
if "smooth_hmm.csv" in dircontent
  fbs_barber = readcsv("smooth_hmm.csv") #read in the ideal answers
else
  fbs_barber = readcsv(joinpath("test","smooth_hmm.csv"))
end
if "smooth_hmm.csv" in dircontent
  filter_barber = readcsv("filter_hmm.csv") # read in the ideal answers
else
  filter_barber = readcsv(joinpath("test","filter_hmm.csv"))
end

# Discrete model
A = [0.5 0.0 0.0;0.3 0.6 0.0;0.2 0.4 1.0] #transition probabilities
B = [0.7 0.4 0.8;0.3 0.6 0.2] # emission probabilities

mod1 = HMM.hmm(A, B) # create the HMM object
initial = [0.9, 0.1, 0.0] # initial state distribution
evidence = [1,1,2,1,2,1,2] # evidence/observations

filter_me = HMM.forward(mod1, initial, evidence)


fbs_me = zeros(length(initial), length(evidence))
for k=1:length(evidence)
  fbs_me[:, k] = HMM.smooth(mod1, initial, evidence, k) # works!
end

vtb_me = HMM.viterbi(mod1, initial, evidence) # works!
vtb_barber = [1, 3, 3, 3, 3, 3, 3] # Barber's answer

pstates, pevidence = HMM.prediction(mod1, initial, evidence) # No test for this - not implemented by barber

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
# Smoother Inference
smooth_handler(r::Test.Success) = println("Successful smoothing test!")
smooth_handler(r::Test.Failure) = error("Failure with the moothing test: $(r.expr)")
smooth_handler(r::Test.Error) = rethrow(r)
Test.with_handler(smooth_handler) do
  @test maximum(abs(fbs_me - fbs_barber)) < 1e-4
end
