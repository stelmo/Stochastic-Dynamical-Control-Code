% Du Reactor Linearisation
clear all
close all
clc
format short
options = optimset('Display','off');

%% Specify the reactor dynamics
phi = 0.072;
q = 1.0;
beta = 8.0;
delta = 0.3;
lambda = 20.0;
x1f = 1.0;
x2f = 0.0;
u = 0.0; % controller input!
k = @(x) exp(x/(1+x/lambda)); 

fx1 = @(x) -phi*k(x(2))-q;
fx2 = @(x) -phi*x(1)*k(x(2))*(lambda/(x(2)+lambda))^2;
gx1 = @(x) beta*phi*k(x(2));
gx2 = @(x) beta*phi*x(1)*k(x(2))*(lambda/(x(2)+lambda))^2-(q+delta);

J = @(x) [fx1(x) fx2(x); gx1(x) gx2(x)]; % the jacobian
sys = @(t, x) [-phi*x(1)*k(x(2))+q*(x1f-x(1)); 
              beta*phi*x(1)*k(x(2))-(q+delta)*x(2)+delta*u+q*x2f]; % For ODE solving
          
guess1 = [0.23, 4.71];
guess2 = [0.55, 2.75];
guess3 = [0.85, 0.88];
[xss1, ~] = fsolve(@(x) sys(1.0, x), guess1, options);
[xss2, ~] = fsolve(@(x) sys(1.0, x), guess2, options);
[xss3, ~] = fsolve(@(x) sys(1.0, x), guess3, options);      

J1 = J(xss1);
% e1 = eig(J1) %stability for critical point
J2 = J(xss2);
% e2 = eig(J2) %stability for critical point
J3 = J(xss3);
% e3 = eig(J3) %stability for critical point

lin1 = @(x, Q) J1*x - J1*xss1' + [0.0; delta]*Q;
lin2 = @(x, Q) J2*x - J2*xss2' + [0.0; delta]*Q;
lin3 = @(x, Q) J3*x - J3*xss3' + [0.0; delta]*Q;

%% Approximation 1
tend = 10;
init1 = [0.15; 4.5];
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
hold off
legend([p1, p1lin],'Nonlinear ODE','Linearised ODE');
ylabel('Concentration [C_A]','fontsize', 18)
set(gca,'fontsize', 18);

subplot(2,1,2)
p2 = plot(t1,y1(:,2),'b');
hold on
p2lin = plot(t1s, y1lin(2,:),'r');
hold off
legend([p2, p2lin],'Nonlinear ODE','Linearised ODE');
xlabel('Time [min]','fontsize', 18)
ylabel('Temperature [K]','fontsize', 18)
set(gca,'fontsize', 18);

%% Approximation 2
tend = 5;
init1 = [0.54; 2.7];
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
hold off
legend([p1, p1lin],'Nonlinear ODE','Linearised ODE');
ylabel('Concentration [C_A]','fontsize', 18)
set(gca,'fontsize', 18);

subplot(2,1,2)
p2 = plot(t1,y1(:,2),'b');
hold on
p2lin = plot(t1s, y1lin(2,:),'r');
hold off
legend([p2, p2lin],'Nonlinear ODE','Linearised ODE');
xlabel('Time [min]','fontsize', 18)
ylabel('Temperature [K]','fontsize', 18)
set(gca,'fontsize', 18);

%% Approximation 3
tend = 5.0;
init1 = [0.99; 0.99];
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
hold off
legend([p1, p1lin],'Nonlinear ODE','Linearised ODE');
ylabel('Concentration [C_A]','fontsize', 18)
set(gca,'fontsize', 18);

subplot(2,1,2)
p2 = plot(t1,y1(:,2),'b');
hold on
p2lin = plot(t1s, y1lin(2,:),'r');
hold off
legend([p2, p2lin],'Nonlinear ODE','Linearised ODE');
xlabel('Time [min]','fontsize', 18)
ylabel('Temperature [K]','fontsize', 18)
set(gca,'fontsize', 18);