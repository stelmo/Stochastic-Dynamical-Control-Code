# Wood Berry Column without dead time.
# The model is taken from: Terminal composition control of a binary distillation
# column by R.K. Wood and M.W. Berry. From Chemical Engineering Science 1973.
# See the folder: literature/control on dropbox
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
b = zeros(4)

sigmaPlant = 0.1
sigmaMeasure = 10.0
Q = sigmaPlant^2*eye(4) # modelling errors
R = sigmaMeasure^2*eye(2) # measurement noise

wb = LLDS_functions.LLDS(A,B,b,C,D,Q,R)

T = 1000; time = [1:1:T] # simulation time
hidden_states = zeros(4, T) # hidden states
visible_states = zeros(2, T) # noisy visible states
visible_states_nn = zeros(2, T) # noiseless visible states

# Controller action
controller_input = zeros(2, T)
controller_input[1, 100:500] = 1.0
controller_input[2, 250:600] = 1.0

# Simulate Model
hidden_states[:,1] = 0.0 # starts at zero
hidden_states[:, 1], visible_states[:, 1], visible_states_nn[:, 1] = LLDS_functions.step(hidden_states[:, 1], controller_input[:, 1], wb)
for t=2:T
  hidden_states[:, t], visible_states[:, t], visible_states_nn[:, t] = LLDS_functions.step(hidden_states[:, t-1], controller_input[:, t-1], wb)
end

# Filter
init_mean = zeros(4) # in deviation therefore I know the starting point is at  0
init_covar = eye(4) #vague prior
filtermeans = zeros(4, T)
filtercovar = zeros(4, 4, T)
filtermeans[:, 1], filtercovar[:, :, 1] = LLDS_functions.init_filter(init_mean, init_covar, controller_input[:,1], visible_states[:,1], wb)
for t=2:T
  filtermeans[:, t], filtercovar[:, :, t] = LLDS_functions.step_filter(filtermeans[:, t-1], filtercovar[:,:, t-1], controller_input[:,t], visible_states[:,t], wb)
end

inferred_visible_states = zeros(2, T) # not sure about this
for t=1:T
  inferred_visible_states[:,t] = wb.C*filtermeans[:, t] + wb.D*controller_input[:, t]
end

# Prediction
pstart = 210 # start predicting here (inclusive)
pend = T # prediction horizon
pred_us = controller_input[:, pstart:pend]
pmeans, pcovars = LLDS_functions.predict_visible(filtermeans[:, pstart-1], filtercovar[:,:, pstart-1], pred_us, wb)


# Plotting
figure(1) # Filtering
vs1_nn, = plot(time, visible_states_nn[1,:]', "r")
vs2_nn, = plot(time, visible_states_nn[2,:]', "b")
vs1, = plot(time, visible_states[1,:]', "rx")
vs2, = plot(time, visible_states[2,:]', "bx")
ivs1, = plot(time[1:10:end], inferred_visible_states[1, 1:10:end]', "r>")
ivs2, = plot(time[1:10:end], inferred_visible_states[2, 1:10:end]', "b>")
pvs1, = plot(time[pstart:10:pend], pmeans[1, 1:10:end]', "rs")
pvs2, =plot(time[pstart:10:pend], pmeans[2, 1:10:end]', "bs")
legend([vs1_nn,vs2_nn, vs1, vs2, ivs1, ivs2, pvs1, pvs2], ["X_D","X_B","Measured X_D","Measured X_B","Filtered X_D","Filtered X_B","Predicted X_D","Predicted X_B"])
