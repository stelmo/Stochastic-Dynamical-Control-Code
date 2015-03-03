% CSTR Linearisation
clear all
close all
clc
format short
options = optimset('Display','off');

%% Specify the reactor parameters
V = 0.1; %m3 
R = 8.314; %kJ/kmol.K
CA0 = 1.0; %kmol/m3
TA0 = 310.0; %K
Q = 0.0; %kJ/min
dH = -4.78e4; %kJ/kmol
k0 = 72.0e9; %1/min
E = 8.314e4; %kJ/kmol
Cp = 0.239; %kJ/kgK
rho = 1000.0; %kg/m3
F = 100e-3; %m3/min

%% Reactor dynamics
% x(1) = CA
% x(2) = TR
sys = @(t, x) [F/V*CA0 - (F/V + k0*exp(-E/(R*x(2))))*x(1);
                F/V*TA0 - F/V*x(2) - dH/(rho*Cp)*k0*exp(-E/(R*x(2)))*x(1)+Q/(rho*Cp*V)];

J = @(x) [-F/V-k0*exp(-E/(R*x(2))),-x(1)*k0*exp(-E/(R*x(2)))*(E/(R*x(2)^2));-dH/(rho*Cp)*k0*exp(-E/(R*x(2))),-(F/V + dH/(rho*Cp)*k0*exp(-E/(R*x(2)))*(E/(R*x(2)^2))*x(1))];

guess1 = [0.0, 509.];
guess2 = [0.57, 395.];
guess3 = [0.99, 310.];
[xss1, ~] = fsolve(@(x) sys(1.0, x), guess1, options);
[xss2, ~] = fsolve(@(x) sys(1.0, x), guess2, options);
[xss3, ~] = fsolve(@(x) sys(1.0, x), guess3, options);      

J1 = J(xss1);
% e1 = eig(J1) $stability for critical point
J2 = J(xss2);
% e2 = eig(J2) $stability for critical point
J3 = J(xss3);
% e3 = eig(J3) $stability for critical point

lin1 = @(x, Q) J1*x - J1*xss1' + [0.0; 1./(rho*Cp*V)]*Q;
lin2 = @(x, Q) J2*x - J2*xss2' + [0.0; 1./(rho*Cp*V)]*Q;
lin3 = @(x, Q) J3*x - J3*xss3' + [0.0; 1./(rho*Cp*V)]*Q;

%% Check Approximations
init1 = [0.57; 390];
tend = 5;
%% Approximation 1
[t1, y1] = ode45(sys, [0 tend], init1);

h1 = 0.001;
t1s = 0:h1:tend;
y1lin = zeros(2,length(t1s));
y1lin(:, 1) = init1;
for k=2:length(t1s)
    y1lin(:,k) = y1lin(:,k-1) + h1*lin1(y1lin(:,k-1),0.0);
end

figure(1)
subplot(2,1,1)
p1 = plot(t1,y1(:,1),'b');
hold on
p1lin = plot(t1s, y1lin(1,:),'r');
ylim([0.0, 1.0])
hold off
legend([p1, p1lin],'Nonlinear ODE','Linearised ODE');
ylabel('Concentration [C_A]','fontsize', 18)
set(gca,'fontsize', 18);

subplot(2,1,2)
p2 = plot(t1,y1(:,2),'b');
hold on
p2lin = plot(t1s, y1lin(2,:),'r');
ylim([250, 550])
hold off
legend([p2, p2lin],'Nonlinear ODE','Linearised ODE');
xlabel('Time [min]','fontsize', 18)
ylabel('Temperature [K]','fontsize', 18)
set(gca,'fontsize', 18);

%% Approximation 2
% init1 = [0.54; 395.];
[t1, y1] = ode45(sys, [0 tend], init1);

h1 = 0.0001;
t1s = 0:h1:tend;
y1lin = zeros(2,length(t1s));
y1lin(:, 1) = init1;
for k=2:length(t1s)
    y1lin(:,k) = y1lin(:,k-1) + h1*lin2(y1lin(:,k-1),0.0);
end

figure(2)
subplot(2,1,1)
p1 = plot(t1,y1(:,1),'b');
hold on
p1lin = plot(t1s, y1lin(1,:),'r');
ylim([0.0, 1.0])
hold off
legend([p1, p1lin],'Nonlinear ODE','Linearised ODE');
ylabel('Concentration [C_A]','fontsize', 18)
set(gca,'fontsize', 18);

subplot(2,1,2)
p2 = plot(t1,y1(:,2),'b');
hold on
p2lin = plot(t1s, y1lin(2,:),'r');
ylim([250, 550])
hold off
legend([p2, p2lin],'Nonlinear ODE','Linearised ODE');
xlabel('Time [min]','fontsize', 18)
ylabel('Temperature [K]','fontsize', 18)
set(gca,'fontsize', 18);

%% Approximation 3
% init1 = [0.65; 400];
[t1, y1] = ode45(sys, [0 tend], init1);

h1 = 0.001;
t1s = 0:h1:tend;
y1lin = zeros(2,length(t1s));
y1lin(:, 1) = init1;
for k=2:length(t1s)
    y1lin(:,k) = y1lin(:,k-1) + h1*lin3(y1lin(:,k-1),0.0);
end

figure(3)
subplot(2,1,1)
p1 = plot(t1,y1(:,1),'b');
hold on
p1lin = plot(t1s, y1lin(1,:),'r');
ylim([0.0, 1.0])
hold off
legend([p1, p1lin],'Nonlinear ODE','Linearised ODE');
ylabel('Concentration [C_A]','fontsize', 18)
set(gca,'fontsize', 18);

subplot(2,1,2)
p2 = plot(t1,y1(:,2),'b');
hold on
p2lin = plot(t1s, y1lin(2,:),'r');
ylim([250, 550])
hold off
legend([p2, p2lin],'Nonlinear ODE','Linearised ODE');
xlabel('Time [min]','fontsize', 18)
ylabel('Temperature [K]','fontsize', 18)
set(gca,'fontsize', 18);
