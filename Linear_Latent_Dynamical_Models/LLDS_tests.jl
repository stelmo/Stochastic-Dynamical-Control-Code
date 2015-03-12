## Linear Latent Dynamic Models
# u   u
# |   |
# h - h -> ...
# |   |
# o   o
# Tests for the Kalman type models.
# All tests are compared to code supplied by:
# Bayesian Reasoning and Machine Learning
# by David Barber
# Website: http://www0.cs.ucl.ac.uk/staff/d.barber/brml/
# The example is taken from Example 24.4
# The functions were all compared to his demo program in Matlab:
# demoLDSTracking.m

using Base.Test
using Distributions

import LLDS_functions

# Specify System
model = begin
  dt = 0.1
  T = 400
  A = eye(6) # state transition
  A[1,5] = dt
  A[2,1] = dt
  A[3,6] = dt
  A[4,3] = dt
  B = zeros(6,1) # input transition
  b = zeros(6) # constant in the transition function
  C = zeros(2,6) # state observation
  C[1,2] = 1.0
  C[2,4] = 1.0
  sigmaQ = 0.00001 # standard deviation of process noise
  sigmaR = 50.0 # standard deviation of measurement noise
  Q = sigmaQ^2*eye(6) # process noise covariance
  R = sigmaR^2*eye(2) # measurement noise covariance

  LLDS_functions.LLDS(A, B, b, C, Q, R)
end

# Specify initial conditions
init_covar = eye(6) # vague prior covar
init_mean = zeros(6) # vague prior mean

# Read in correct data
dircontent = readdir()
if "visiblestates.csv" in dircontent
  visiblestates = readcsv("visiblestates.csv") #read in the ideal answers
else
  visiblestates = readcsv(string(pwd(),"/Linear_Latent_Dynamical_Models/visiblestates.csv"))
end
if "filtercovar.csv" in dircontent
  filtercovar = reshape(readcsv("filtercovar.csv"), 6,6,T) #read in the ideal answers
else
  filtercovar = reshape(readcsv(string(pwd(),"/Linear_Latent_Dynamical_Models/filtercovar.csv")), 6, 6, T)
end
if "filtermeans.csv" in dircontent
  filtermeans = readcsv("filtermeans.csv") #read in the ideal answers
else
  filtermeans = readcsv(string(pwd(),"/Linear_Latent_Dynamical_Models/filtermeans.csv"))
end
if "smoothcovar.csv" in dircontent
  smoothedcovar = reshape(readcsv("smoothcovar.csv"), 6,6, T) #read in the ideal answers
else
  smoothedcovar = reshape(readcsv(string(pwd(),"/Linear_Latent_Dynamical_Models/smoothcovar.csv")),6,6,T)
end
if "smoothmeans.csv" in dircontent
  smoothedmeans = readcsv("smoothmeans.csv") #read in the ideal answers
else
  smoothedmeans = readcsv(string(pwd(),"/Linear_Latent_Dynamical_Models/smoothmeans.csv"))
end

# Filter
ucontrol = zeros(1) # no control so this is really only a dummy variable.
filtermeans_own = zeros(6, T)
filtercovar_own = zeros(6, 6, T)
filtermeans_own[:, 1], filtercovar_own[:, :, 1] = LLDS_functions.init_filter(init_mean, init_covar, visiblestates[:,1], model)
for t=2:T
  filtermeans_own[:, t], filtercovar_own[:, :, t] = LLDS_functions.step_filter(filtermeans_own[:, t-1], filtercovar_own[:,:, t-1], ucontrol, visiblestates[:,t], model)
end

# Smoothed
ucontrols = zeros(1, T) # no control so this is really only a dummy variable.
smoothedmeans_own = zeros(6, T)
smoothedcovar_own = zeros(6, 6, T)

smoothedmeans_own, smoothedcovar_own = LLDS_functions.smooth(filtermeans_own, filtercovar_own, ucontrols, model)

# Run the tests
tol = 0.01
function unroll(a::Array{Float64, 3}, b::Array{Float64, 3}, T::Int64)
  # Finds maximum difference between two 3d matrices.
  maxdiff = 0.0
  for t=1:T
    tempdiff = maximum(abs(a[:,:,t] - b[:,:,t]))
    (maxdiff < tempdiff) && (maxdiff = tempdiff)
  end
  return maxdiff
end

# Filter Inference
filtermean_handler(r::Test.Success) = println("Successful filter mean test for LLDS!")
filtermean_handler(r::Test.Failure) = error("Failure with the filter mean test  for LLDS: $(r.expr)")
filtermean_handler(r::Test.Error) = rethrow(r)
Test.with_handler(filtermean_handler) do
  @test maximum(abs(filtermeans_own - filtermeans)) < tol
end
filtercovar_handler(r::Test.Success) = println("Successful filter covariance test for LLDS!")
filtercovar_handler(r::Test.Failure) = error("Failure with the filter covariance test for LLDS: $(r.expr)")
filtercovar_handler(r::Test.Error) = rethrow(r)
Test.with_handler(filtercovar_handler) do
  @test unroll(filtercovar_own, filtercovar, T) < tol
end

# Smoothing Inference
smoothingmean_handler(r::Test.Success) = println("Successful smoothing mean test for LLDS!")
smoothingmean_handler(r::Test.Failure) = error("Failure with the smoothing mean test for LLDS: $(r.expr)")
smoothingmean_handler(r::Test.Error) = rethrow(r)
Test.with_handler(smoothingmean_handler) do
  @test maximum(abs(smoothedmeans_own - smoothedmeans)) < tol
end
smoothingcovar_handler(r::Test.Success) = println("Successful smoothing covariance test for LLDS!")
smoothingcovar_handler(r::Test.Failure) = error("Failure with the smoothing covariance test for LLDS: $(r.expr)")
smoothingcovar_handler(r::Test.Error) = rethrow(r)
Test.with_handler(smoothingcovar_handler) do
  @test unroll(smoothedcovar_own, smoothedcovar, T) < tol
end
