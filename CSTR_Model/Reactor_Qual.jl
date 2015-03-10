# Qualitative Analysis of the CSTR

using PyPlot
import Reactor_functions
reload("Reactor_functions.jl")

# First, lets inspect the steady state values for different parameters
cstr1 = begin
  V = 0.1 #m3
  R = 8.314 #kJ/kmol.K
  CA0 = 1.0 #kmol/m3
  TA0 = 310.0 #K
  dH = -4.78e4 #kJ/kmol
  k0 = 72.0e9 #1/min
  E = 8.314e4 #kJ/kmol
  Cp = 0.239 #kJ/kgK
  rho = 1000.0 #kg/m3
  F = 100e-3 #m3/min
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

cstr2 = begin
  V = 0.1 #m3
  R = 8.314 #kJ/kmol.K
  CA0 = 1.0 #kmol/m3
  TA0 = 310.0 #K
  dH = -4.78e4 #kJ/kmol
  k0 = 72.0e8 #1/min
  E = 8.314e4 #kJ/kmol
  Cp = 0.239 #kJ/kgK
  rho = 1000.0 #kg/m3
  F = 100e-3 #m3/min
  Reactor_functions.Reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

N = 100
Ts = linspace(250, 550, N) # temperature range
qrs1 = zeros(N) # heat removals
qrs2 = zeros(N) # heat removals
qrs3 = zeros(N) # heat removals
qgs1 = zeros(N) # heat generations
qgs2 = zeros(N) # heat generations

for k=1:N
  qrs1[k] = Reactor_functions.QR(Ts[k], -1250.0, cstr1)
  qrs2[k] = Reactor_functions.QR(Ts[k], 0.0, cstr1)
  qrs3[k] = Reactor_functions.QR(Ts[k], 950.0, cstr1)
  qgs1[k] = Reactor_functions.QG(Ts[k], cstr1)
  qgs2[k] = Reactor_functions.QG(Ts[k], cstr2)
end

figure(1)
q1, = plot(Ts, qrs1, "b", linewidth=3)
q2, = plot(Ts, qrs2, "g", linewidth=3)
q3, = plot(Ts, qrs3, "r", linewidth=3)
plot(Ts, qgs1, "k", linewidth=3)
# plot(Ts, qgs2, "k--", linewidth=3)
legend([q1,q2,q3],["Q=-1250 kJ/min","Q=0 kJ/min","Q=950 kJ/min"], loc="best")
ylim([0, maximum([qrs, qgs])])
xlabel("Steady State Temperature [K]")
ylabel("Heat Removal [K/min]")
