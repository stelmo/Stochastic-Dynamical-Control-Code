# Simulation of the simplified Tennessee Eastman Challenge Process

include("Plant.jl")
using PyPlot
ioff()

ste = Plant.plant()

N = 50000 # time discretisation
ts = linspace(0, 60, N)
h = ts[2]-ts[1]
delay = round(0.1/h) # steps required to move 6 min

# Storage matrices
xs = zeros(8, N) # states
xs[:,1] = [46.2281, 19.2127, 32.9168, 108.436, 60.3487, 25.2741, 27.9094, 46.1511]
us = zeros(4, N) # inputs
us[1, :] = 60.3487 # u1 initial
us[2, :] = 25.2741 # u2 initial
us[3, :] = 27.9094 # u3 initial
us[4, :] = 44.1797 # u4 initial
ys = zeros(10, N) # outputs
ys[:,1] = Plant.getOutputs(xs[:,1])
es = zeros(4, N) # signal errors

# Controller specifications
F4_sp = 100.
P_sp = 2800.
yA3_sp = 47.
es[1, 1] = F4_sp - ys[4, 1] # Connects F4 to u1
es[2, 1] = yA3_sp - ys[7, 1] # Connects ya3 to u2
es[3, 1] = P_sp - ys[5, 1] # Connects P to u3
P_max = 2900
es[4, 1] = P_max - ys[5, 1]
F4_adj = 0.0
for k=2:N

   # Tests
   if k > N/5.
      yA3_sp = 60.
   end

   # Model Update
   xs[:, k] = Plant.stepModelRK(xs[:,k-1], us[:,k-1], h, ste)

   # Measurement Update
   ys[:, k] = Plant.getOutputs(xs[:,k])
   # Remember pure time delays for composition measurement
   if float(k % delay) != 0.0
      # Measurement delay from the chromatograph
      ys[7, k] = ys[7, k-1]
      ys[8, k] = ys[8, k-1]
      ys[9, k] = ys[9, k-1]
   end

   # Override
   es[4, k] = P_max - ys[5, k]
   F4_adj = Plant.override(F4_adj, es[4, k], es[4, k-1], 0.7, 3.0)
   F4_sp_in = F4_sp + F4_adj

   # Error update
   es[1, k] = F4_sp_in - ys[4, k] # Connects F4 to u1
   es[2, k] = yA3_sp - ys[7, k] # Connects ya3 to u2
   es[3, k] = P_sp - ys[5, k] # Connects P to u3

   # Control System
   if float(k % delay) == 0.0
      # Discrete controllers
      us[1, k] = Plant.controller(us[1, k-1], es[1, k], es[1, k-1], 0.1, 1.0)
      us[2, k] = Plant.controller(us[2, k-1], es[2, k], es[2, k-1], 2.0, 3.0)
      us[3, k] = Plant.controller(us[3, k-1], es[3, k], es[3, k-1], -0.5, 1.5)
   else
      us[1, k] = us[1, k-1]
      us[2, k] = us[2, k-1]
      us[3, k] = us[3, k-1]
   end

end

# Plant.plotStates(ts, xs)
Plant.plotInputs(ts, us)
Plant.plotOutputs(ts, ys)

plt.show()
