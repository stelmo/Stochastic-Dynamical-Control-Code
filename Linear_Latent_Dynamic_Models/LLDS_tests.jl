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

using PyPlot

import LLDS_functions
reload("LLDS_functions.jl")

# Specify System
dt = 0.1
T = 400
A = eye(6) # state transition
A[1,5] = dt
A[2,1] = dt
A[3,6] = dt
A[4,3] = dt
B = zeros(6,1) # input transition
C = zeros(2,6) # state observation
C[1,2] = 1.0
C[2,4] = 1.0
D = zeros(2,1) # input observation

sigmaQ = 0.00001 # standard deviation of process noise
sigmaR = 50.0 # standard deviation of measurement noise
Q = sigmaQ^2*eye(6) # process noise covariance
R = sigmaR^2*eye(2) # measurement noise covariance

model = LLDS_functions.LLDS{Array{Float64,2}}(A, B, C, D, Q, R)

# Specify initial conditions
init_covar = eye(6) # vague prior covar
init_mean = zeros(6) # vague prior mean

# Read in correct data
hiddenstates = readcsv("hiddenstates.csv")
visiblestates = readcsv("visiblestates.csv")

filtercovar = readcsv("filtercovar.csv")
filtermeans = readcsv("filtermeans.csv")

smoothcovar = readcsv("smoothcovar.csv")
smoothmeans = readcsv("smoothmeans.csv")


# Plotting
figure(1)
plot(hiddenstates[2,1:10:end], hiddenstates[4,1:10:end], "ro")
plot(visiblestates[1,:], visiblestates[2,:], "g.")
plot(filtermeans[2, 1:10:end], filtermeans[4,1:10:end], "c>")
plot(smoothmeans[2, 1:10:end], smoothmeans[4,1:10:end], "bs")
