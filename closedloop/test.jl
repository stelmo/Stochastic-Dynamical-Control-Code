# temp test file
include("../init_sys.jl")

using LLDS
using Reactor
using Ellipse
using LQR

# Specify the nonlinear model
cstr = begin
  V = 5.0 #m3
  R = 8.314 #kJ/kmol.K
  CA0 = 1.0 #kmol/m3
  TA0 = 310.0 #K
  dH = -4.78e4 #kJ/kmol
  k0 = 72.0e7 #1/min
  E = 8.314e4 #kJ/kmol
  Cp = 0.239 #kJ/kgK
  rho = 1000.0 #kg/m3
  F = 100e-3 #m3/min
  Reactor.reactor(V, R, CA0, TA0, dH, k0, E, Cp, rho, F)
end

init_state = [0.50; 400]
h = 0.1 # time discretisation
tend = 150.0 # end simulation time
ts = [0.0:h:tend]
N = length(ts)
xs = zeros(2, N)
xs[:,1] = init_state
ys = zeros(2, N) # only one measurement

xspace = [0.0, 1.0]
yspace = [250, 650]


# Specify the linear model
linsystems = Reactor.getLinearSystems(nX, nY, xspace, yspace, h, cstr_model)
QQ = zeros(2, 2)
QQ[1] = 1.0
RR = 1.0
H = [1. 0.]

NN = length(linsystems)
KK =
for n=1:NN

end
A = linsystems[2].A
B = linsystems[2].B
b = linsystems[2].b
C = eye(2)

DARE = LQR.dare(A, B, QQ, RR)
K = LQR.lqr(A, B, QQ, RR)

xoff, uoff = LQR.offset(A, B, C, H, ysp)
