# Wood Berry Column without dead time.
# The model is taken from: Terminal composition control of a binary distillation
# column by R.K. Wood and M.W. Berry. From Chemical Engineering Science 1973.
# The model is:
# [xD] = [P11 P12] [R]
# [xB] = [P21 P22] [S]
# with P11 = 12.8/(16.7s+1) and P12 = -18.9/(21.0*s+1)
# and P21 = 6.6/(10.9*s+1) and P22 = -19.4/(14.4*s+1)

using PyPlot
using LLDS_functions

A = readcsv("WB_A.csv")
B = readcsv("WB_B.csv")
C = readcsv("WB_C.csv")
D = readcsv("WB_D.csv")

sigmaPlant = 0.01
sigmaMeasure = 1.0
Q = sigmaPlant^2*eye(4) # modelling errors
R = sigmaMeasure^2*eye(2) # measurement noise

wb = LLDS_functions.LLDS(A,B,C,D,Q,R)

T = 1000; time = [1:1:T] # simulation time
hidden_states = zeros(4, T) # hidden states
visible_states = zeros(2, T) # noisy visible states
visible_states_nn = zeros(2, T) # noiseless visible states

controller_input = zeros(2, T)
controller_input[1, 100:500] = 1.0
controller_input[2, 250:600] = 1.0

for t=2:T
  hidden_states[:, t], visible_states[:, t], visible_states_nn[:, t] = LLDS_functions.step(hidden_states[:, t-1], controller_input[:, t-1], wb)
end

mean_states = zeros(4, T)
covar_states = sigmaPlant^2*eye(4, 4*T)

# Plotting
figure(1)
plot(time, visible_states_nn[1,:]', "r")
plot(time, visible_states_nn[2,:]', "b")
plot(time, visible_states[1,:]', "rx")
plot(time, visible_states[2,:]', "bx")
