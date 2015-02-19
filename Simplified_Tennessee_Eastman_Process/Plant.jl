module Plant
# A model for the simplified Tennessee Eastman Challenge problem

using DataFrames
using Distributions
using PyPlot

type plant
   R :: Float64 # Ideal gas constant
   T :: Float64 # Temperature
   V :: Float64 # Volume
   rhoL :: Float64 # Liquid density of D
   k0 :: Float64 # Reaction constant
   F1max :: Float64 # Feed 1 max molar flow rate
   F2max :: Float64 # Feed 2 max molar flow rate
   Cv3 :: Float64 # Valve constant
   Cv4 :: Float64 # Valve constant
   tauV :: Float64 # Valve time constant
   x4bar :: Float64 # SS of molar hold up D
   Kc :: Float64 # Proportional constant
   yA1 :: Float64 # Mole frac in Feed 1 of A
   yB1 :: Float64 # Mole frac in Feed 1 of B
   yC1 :: Float64 # Mole frac in Feed 1 of A

   function plant() # create the ideal plant
   new(8.314, 373.0, 122.0, 8.3, 0.00117, 330.46, 22.46, 0.00352, 0.0417, 0.002777, 47.0302, -1.4, 0.485, 0.005, 0.510)
   end
end

function reaction(nA, nB, nC, nD, model)
   # Reaction term
   return  model.k0 * ( ((nA*model.R*model.T)/(model.V - (nD/model.rhoL)) )^1.2 ) * (((nC*model.R*model.T)/(model.V - (nD/model.rhoL)) )^0.4)
end

function valveType2(nA, nB, nC, nD, xI, CvI, model)
   # General valve term of the second type. NB: excludes the species specific term!
   return CvI * sqrt(((nA+nB+nC)*model.R*model.T) / (model.V-(nD/model.rhoL)) - 100.0) * xI # modified
end

function sysfunc(y, u, model)
   # The Forward Euler system of nonlinear ODEs

   fy = zeros(8)

   # variables order
   nA = 1
   nB = 2
   nC = 3
   nD = 4
   x1 = 5
   x2 = 6
   x3 = 7
   x4 = 8

   #x1 = 60.95327
   #x2 = 25.022322
   #x3 = 39.25777
   #x4 = 47.03024

   # Function - mistakes in the paper: F3, F4 need to me multiplied by a factor of 100
   fy[1] = model.yA1*model.F1max*(y[x1]/100.0) + model.F2max*(y[x2]/100.0) - (y[nA]/(y[nA]+y[nB]+y[nC]))* valveType2(y[nA], y[nB], y[nC], y[nD], y[x3], model.Cv3, model) - reaction(y[nA], y[nB], y[nC], y[nD], model)
   fy[2] = model.yB1*model.F1max*(y[x1]/100.0) - (y[nB]/(y[nA]+y[nB]+y[nC]))* valveType2(y[nA], y[nB], y[nC], y[nD], y[x3], model.Cv3, model)
   fy[3] = model.yC1*model.F1max*(y[x1]/100.0) - (y[nC]/(y[nA]+y[nB]+y[nC]))* valveType2(y[nA], y[nB], y[nC], y[nD], y[x3], model.Cv3, model) - reaction(y[nA], y[nB], y[nC], y[nD], model)
   fy[4] = reaction(y[nA], y[nB], y[nC], y[nD], model) - valveType2(y[nA], y[nB], y[nC], y[nD], y[x4], model.Cv4, model)
   fy[5] = (1/model.tauV)*(u[1] - y[x1])
   fy[6] = (1/model.tauV)*(u[2] - y[x2])
   fy[7] = (1/model.tauV)*(u[3] - y[x3])
   fy[8] = (1/model.tauV)*( (model.x4bar + model.Kc*(u[4] - (y[nD]/model.rhoL)/0.30)) - y[x4] )

   return fy
end

function stepModelFE(yold, uold, h, model)
   # Forward Euler integration step

   ynew = yold + h*sysfunc(yold, uold, model)
   return ynew
end

function stepModelRK(yold, uold, h, model)
   # Runge Kutta solution
   # There is no explicit time dependency. It could be
   # argued that the forcing functions are time dependent
   # but this would be strange to apply discretely? Look into it!
   k1 = sysfunc(yold, uold, model)
   k2 = sysfunc(yold + h*k1*0.5, uold, model)
   k3 = sysfunc(yold + h*k2*0.5, uold, model)
   k4 = sysfunc(yold + h*k3, uold, model)
   ynew = yold + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
   return ynew
end

function addNoise(ynow)
   # Adds random noise to each channel.
   C = 0.001*eye(8)
   C[2,2] = 0.001
   ymod = ynow + rand(MvNormal(C))
   return ymod
end

function createDF(ts, ys)
   # Returns a data frame for plotting.
   n = length(ts)
   df = DataFrame() # row#, time, data, key
   k1 = fill("NA", n)
   k2 = fill("NB", n)
   k3 = fill("NC", n)
   k4 = fill("ND", n)
   k5 = fill("X1", n)
   k6 = fill("X2", n)
   k7 = fill("X3", n)
   k8 = fill("X4", n)

   df[:time] = [ts, ts, ts, ts, ts, ts, ts, ts]
   df[:data] = [ys[1,:][:], ys[2,:][:], ys[3,:][:], ys[4,:][:], ys[5,:][:], ys[6,:][:], ys[7,:][:], ys[8,:][:]]
   df[:key] = [k1, k2, k3, k4, k5, k6, k7, k8]

   return df
end

function plotInputs(ts, us)

  #us = round(us, 2)
  figure()
  suptitle("Inputs")

  subplot(141)
  plot(ts, us[1, :][:])
  title("Valve Signal: \$X_1\$")
  xlabel("Time [hours]")
  ylim(0., 100.)

  subplot(142)
  plot(ts, us[2, :][:])
  title("Valve Signal: \$X_2\$")
  xlabel("Time [hours]")
  ylim(0., 100.)

  subplot(143)
  plot(ts, us[3, :][:])
  title("Valve Signal: \$X_3\$")
  xlabel("Time [hours]")
  ylim(0., 100.)

  subplot(144)
  plot(ts, us[4, :][:])
  title("Set Point \$X_4\$")
  xlabel("Time [hours]")
  ylim(0., 100.)
end

function plotStates(ts, xs)

  #ys = round(ys, 2)
  figure()
  suptitle("States")

  subplot(241)
  x1 = plot(ts, xs[1, :][:])
  title("\$N_A\$")
  ylim(0., 80.)

  subplot(242)
  x2 = plot(ts, xs[2, :][:])
  title("\$N_B\$")
  ylim(0., 30.)

  subplot(243)
  x3 = plot(ts, xs[3, :][:])
  title("\$N_C\$")
  ylim(0., 70.)

  subplot(244)
  x4 = plot(ts, xs[4, :][:])
  title("\$N_D\$")
  ylim(0., 220.)

  subplot(245)
  x5 = plot(ts, xs[5, :][:])
  title("\$X_1\$")
  xlabel("Time [hours]")
  ylim(0., 100.)

  subplot(246)
  x6 = plot(ts, xs[6, :][:])
  title("\$X_2\$")
  xlabel("Time [hours]")
  ylim(0., 100.)

  subplot(247)
  x7 = plot(ts, xs[7, :][:])
  title("\$X_3\$")
  xlabel("Time [hours]")
  ylim(0., 100.)

  subplot(248)
  x8 = plot(ts, xs[8, :][:])
  title("\$X_4\$")
  xlabel("Time [hours]")
  ylim(0., 100.)

end

function plotOutputs(ts, ys)

  figure()
  suptitle("Outputs")

  subplot(251)
  plot(ts, ys[1, :][:])
  title("\$F_1\$")
  ylim(0., 330.46)

  subplot(252)
  plot(ts, ys[2, :][:])
  title("\$F_2\$")
  ylim(0., 22.46)

  subplot(253)
  plot(ts, ys[3, :][:])
  title("\$F_3\$")
  ylim(0., 40.)

  subplot(254)
  plot(ts, ys[4, :][:])
  title("\$F_4\$")
  ylim(0., 200.)

  subplot(255)
  plot(ts, ys[5, :][:])
  title("Pressure")
  ylim(0., 3300.)

  subplot(256)
  plot(ts, ys[6, :][:])
  title("\$V_L\$")
  ylim(0., 100.)

  subplot(257)
  plot(ts, ys[7, :][:])
  title("\$y_{A3}\$")
  ylim(0., 100.)

  subplot(258)
  plot(ts, ys[8, :][:])
  title("\$y_{B3}\$")
  ylim(0., 100.)

  subplot(259)
  plot(ts, ys[9, :][:])
  title("\$y_{C3}\$")
  ylim(0., 100.)

  subplot(2,5,10)
  plot(ts, ys[10, :][:])
  title("Instantaneous Cost")
  ylim(0., 1.)
end

function getOutputs(xs)

   ys = zeros(10)
   ys[1] = 330.46*xs[5]/100.
   ys[2] = 22.46*xs[6]/100.
   # edit valve
   ys[3] = (xs[7])*0.00352*sqrt((xs[1]+xs[2]+xs[3])*8.314*373/(122-xs[4]/8.3)-100.)
   ys[4] = (xs[8])*0.0417*sqrt((xs[1]+xs[2]+xs[3])*8.314*373/(122-xs[4]/8.3)-100.)
   ys[5] = (xs[1]+xs[2]+xs[3])*8.314*373/(122-xs[4]/8.3)
   ys[6] = (xs[4]/8.3)/30*100.
   ys[7] = xs[1]/(xs[1]+xs[2]+xs[3])*100.
   ys[8] = xs[2]/(xs[1]+xs[2]+xs[3])*100.
   ys[9] = xs[3]/(xs[1]+xs[2]+xs[3])*100.
   ys[10] = (ys[3]/ys[4])*(2.206*ys[7]/100. + 6.177*ys[9]/100.)

   return ys
end

function controller(upast, enow, epast, kc, tau)
   # A general discrete PI controller

   unow = upast + kc*(enow - epast + (0.1/tau)*enow)
   #enforce constraints
   if unow > 100.
      unow = 100.
   elseif unow < 0.
      unow = 0.
   end
   return unow
end

function override(upast, enow, epast, kc, tau)
   unow = upast + kc*(enow - epast + (0.1/tau)*enow)
   #enforce constraints
   if unow > 0.
      unow = 0. # saturates at an upperbound of 0
   end
   return unow
end

end # module
