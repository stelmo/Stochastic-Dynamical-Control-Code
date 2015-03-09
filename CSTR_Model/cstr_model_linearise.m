% CSTR Linearisation
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
E = 8.314e4; %kJ/kmol
Cp = 0.239; %kJ/kgK
rho = 1000.0; %kg/m3
F = 100e-3; %m3/min

%% Reactor dynamics
% x(1) = CA
% x(2) = TR
sys = @(t, x, Q) [F/V*CA0 - (F/V + k0*exp(-E/(R*x(2))))*x(1);
                F/V*TA0 - F/V*x(2) - dH/(rho*Cp)*k0*exp(-E/(R*x(2)))*x(1)+Q/(rho*Cp*V)];

J = @(x) [-F/V-k0*exp(-E/(R*x(2))),-x(1)*k0*exp(-E/(R*x(2)))*(E/(R*x(2)^2));-dH/(rho*Cp)*k0*exp(-E/(R*x(2))),-(F/V + dH/(rho*Cp)*k0*exp(-E/(R*x(2)))*(E/(R*x(2)^2))*x(1))];

guess1 = [0.0, 509.];
guess2 = [0.57, 395.];
guess3 = [0.99, 310.];
[xss1, ~] = fsolve(@(x) sys(1.0, x, 0.0), guess1, options);
[xss2, ~] = fsolve(@(x) sys(1.0, x, 0.0), guess2, options);
[xss3, ~] = fsolve(@(x) sys(1.0, x, 0.0), guess3, options);      

J1 = J(xss1);
e1 = eig(J1); %stability for critical point
J2 = J(xss2);
e2 = eig(J2); %stability for critical point
J3 = J(xss3);
e3 = eig(J3); %stability for critical point

dlmwrite('J_ss3.csv', J3, 'precision','%.10f') % write the jacobian
dlmwrite('ss3.csv', xss3, 'precision','%.10f') % write the steady state 

lin1 = @(x, Q) J1*x - J1*xss1' + [0.0; 1./(rho*Cp*V)]*Q;
lin2 = @(x, Q) J2*x - J2*xss2' + [0.0; 1./(rho*Cp*V)]*Q;
lin3 = @(x, Q) J3*x - J3*xss3' + [0.0; 1./(rho*Cp*V)]*Q;

%% Check Approximations
init = [0.5; 410];
tend = 1.0;
h = 0.001;
ts = 0:h:tend;
Qin = 0.0; %input
%% Approximation 1
[tnl, ynl] = ode45(@(t,x) sys(t,x, Qin), [0 tend], init);

y1lin = zeros(2,length(ts));
y1lin(:, 1) = init;
y2lin = zeros(2,length(ts));
y2lin(:, 1) = init;
y3lin = zeros(2,length(ts));
y3lin(:, 1) = init;
ynlfe = zeros(2,length(ts));
ynlfe(:, 1) = init;
for k=2:length(ts)
    y1lin(:,k) = y1lin(:,k-1) + h*lin1(y1lin(:,k-1), Qin);
    y2lin(:,k) = y2lin(:,k-1) + h*lin2(y2lin(:,k-1), Qin);
    y3lin(:,k) = y3lin(:,k-1) + h*lin3(y3lin(:,k-1), Qin);
    ynlfe(:,k) = ynlfe(:,k-1) + h*sys(0.0, ynlfe(:,k-1), Qin); 
end

%% Norm distance
n1 = zeros(1,length(ts));
n2 = zeros(1,length(ts));
err1 = zeros(1,length(ts));
err2 = zeros(1,length(ts));
for k=1:length(ts)
    n1(1,k) = norm(y1lin(:,k)'-xss1);
    n2(1,k) = norm(y2lin(:,k)'-xss2);
    err1(1,k) = norm(y1lin(:,k)'-ynlfe(:,k)');
    err2(1,k) = norm(y2lin(:,k)'-ynlfe(:,k)');
end

figure(1)
subplot(2,1,1)
p1 = plot(tnl,ynl(:,1),'k','LineWidth',3);
hold on
p1lin1 = plot(ts, y1lin(1,:),'r--','LineWidth',3);
p2lin1 = plot(ts, y2lin(1,:),'g--','LineWidth',3);
p3lin1 = plot(ts, y3lin(1,:),'b--','LineWidth',3);
ylim([0.0, 1.0])
hold off
legend([p1, p1lin1],'Nonlinear ODE','Linearised ODE (C_A, T_R)^0_1','Location', 'northeast','Orientation','horizontal');
legend('boxoff')
ylabel('Concentration [kmol.m^{-3}]','fontsize', 18)
set(gca,'fontsize', 18);

subplot(2,1,2)
p2 = plot(tnl,ynl(:,2),'k','LineWidth',3);
hold on
p1lin2 = plot(ts, y1lin(2,:),'r--','LineWidth',3);
p2lin2 = plot(ts, y2lin(2,:),'g--','LineWidth',3);
p3lin2 = plot(ts, y3lin(2,:),'b--','LineWidth',3);
ylim([250, 550])
hold off
legend([p2lin2, p3lin2],'Linearised ODE (C_A, T_R)^0_2','Linearised ODE (C_A, T_R)^0_3','Location', 'southeast','Orientation','horizontal');
legend('boxoff')
xlabel('Time [min]','fontsize', 18)
ylabel('Temperature [K]','fontsize', 18)
set(gca,'fontsize', 18);

figure(2)
subplot(2,1,1)
pedist1 = plot(ts, n1(1,:),'r-','LineWidth',3);
hold on
pedist2 = plot(ts, n2(1,:),'g-','LineWidth',3);
hold off
ylim([0, 200])
ylabel('Euclidean Distance','fontsize', 18)
set(gca,'fontsize', 18);
legend([pedist1, pedist2],'Linearised ODE (C_A, T_R)^0_1', 'Linearised ODE (C_A, T_R)^0_2','Location', 'southeast');

subplot(2,1,2)
perr1 = plot(ts, err1(1,:),'r-','LineWidth',3);
hold on
perr2 = plot(ts, err2(1,:),'g-','LineWidth',3);
hold off
ylim([0, 200])
legend([perr1, perr2],'Linearised ODE (C_A, T_R)^0_1', 'Linearised ODE (C_A, T_R)^0_2','Location', 'northwest');
xlabel('Time [min]','fontsize', 18)
ylabel('Euclidean Error','fontsize', 18)
set(gca,'fontsize', 18);
% 
% h=gcf;
% set(h,'PaperOrientation','landscape');
% set(h,'PaperUnits','normalized');
% set(h,'PaperPosition', [0 0 1 1]);
% print(gcf, '-dpdf', 'cstr_op5.pdf');

