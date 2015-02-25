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

import LLDS_functions

dt = 0.1
T = 400
A = eye(6)
A[1,5] = dt
A[2,1] = dt
A[3,6] = dt
A[4,3] = dt
B = zeros(6,1)
C = zeros(2,6)
C[1,2] = 1.0
C[2,4] = 1.0
D = zeros(2,1)

Q = 1.0e-5*eye(6)
Q[5,5] = 1.0e-2
Q[6,6] = 1.0e-2
R = 50.0*eye(2)

model = LLDS_functions.LLDS{Array{Float64,2}}(A, B, C, D, Q, R)
