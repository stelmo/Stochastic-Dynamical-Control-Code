# Reactor tests: run my implementation of RK on the nonlinear model and compare
# the solution to Matlab's solution.

import Reactor
using Base.Test

dircontent = readdir()
if "state_solutions.csv" in dircontent
  state_solutions = readcsv("state_solutions.csv") #read in the ideal answers
else
  state_solutions = readcsv(joinpath("test","state_solutions.csv"))
end

cstr = begin
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
  Reactor.reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

h = 0.001 # time discretisation
tend = 5. # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
initial_states = [0.57, 395]

xs[:,1] = initial_states
# Loop through the rest of time
for t=2:N
  xs[:, t] = Reactor.run_reactor(xs[:, t-1], 0.0, h, cstr) # actual plant
end

# Run the tests
tol = 0.05
reactor_handler(r::Test.Success) = println("Successful Reactor Integration test!")
reactor_handler(r::Test.Failure) = error("Failure with the Reactor Integration test: $(r.expr)")
reactor_handler(r::Test.Error) = rethrow(r)
Test.with_handler(reactor_handler) do
  @test maximum(abs(state_solutions'-xs)) < tol
end
