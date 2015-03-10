clear all
close all
clc
format long
options = optimset('Display','off');

%% Specify the reactor parameters
V = 0.1; %m3 
R = 8.314; %kJ/kmol.K
CA0 = 1.0; %kmol/m3
TA0 = 310.0; %K
dH = -4.78e4; %kJ/kmol
k0 = 72.0e9; %1/min
E = 8.314e4*2; %kJ/kmol
Cp = 0.239; %kJ/kgK
rho = 1000.0; %kg/m3
F = 100e-3; %m3/min

%% Reactor dynamics
% x(1) = CA
% x(2) = TR
Q = 0.0;
sys = @(t, x) [F/V*CA0 - (F/V + k0*exp(-E/(R*x(2))))*x(1);
                F/V*TA0 - F/V*x(2) - dH/(rho*Cp)*k0*exp(-E/(R*x(2)))*x(1)+Q/(rho*Cp*V)];

guess1 = [0.0, 509.];
[xss1, ~] = fsolve(@(x) sys(1.0, x), guess1, options);

[t, y] = ode45(sys, [0 4], [0.2, 550]);

subplot(2,1,1)
plot(t, y(:,1))
subplot(2,1,2)
plot(t, y(:,2))