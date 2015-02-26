% Du reactor linearisation
clear all
close all
clc
format long
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

F = @(x) [-phi*x(1)*k(x(2))+q*(x1f-x(1)); 
          beta*phi*x(1)*k(x(2))-(q+delta)*x(2)+delta*u+q*x2f]; % For SS solving
      
%% Solve for the steady states
ss1 = [0.85; 0.88];
ss2 = [0.55; 2.75];
ss3 = [0.23; 4.71];
[xss1, fval] = fsolve(F,ss1);
[xss2, fval] = fsolve(F,ss2);
[xss3, fval] = fsolve(F,ss3);
J = @(x) [fx1(x) fx2(x); gx1(x) gx2(x)]; % the jacobian

J1 = J(xss1);
J2 = J(xss2);
J3 = J(xss3);

e1 = eig(J1);
e2 = eig(J2);
e3 = eig(J3);

%% Solve non-linear system
f = @(t, x) [-phi*x(1)*k(x(2))+q*(x1f-x(1)); 
             beta*phi*x(1)*k(x(2))-(q+delta)*x(2)+delta*u+q*x2f]; % For ODE solving

[tout, xout] = ode45(f, [0; 10], [0.8; 0.83]);
%% Solve linear system
Jlin = J1; % set which system to solve
xss = xss1; % set the steady state point
% flin = @(t, x) [Jlin(1,1)*x(1) + Jlin(1,2)*x(2);Jlin(2,1)*x(1) + Jlin(2,2)*x(2)];
% [toutlin, xoutlin] = ode45(flin, [0; 10], [0.8-xss1(1);0.83-xss1(2)]);
% xoutlin = xoutlin + repmat(xss1, 1, length(toutlin))';

flin = @(t, x) [Jlin(1,1)*x(1) + Jlin(1,2)*x(2); 
                Jlin(2,1)*x(1) + Jlin(2,2)*x(2)] - Jlin*xss;
[toutlin, xoutlin] = ode45(flin, [0; 10], [0.8;0.83]);

% test if forward euler is good enough
time = linspace(0.0, 10, 150);
fex = zeros(150,2);
h = time(2);
fex(1, :) = [0.8, 0.83]; 
for k=2:150
   fex(k, :) = fex(k-1,:) + h.*f(0,fex(k-1,:))'; 
end
%% Plot
plot(tout, xout, '--')
hold on
plot(time, fex, '.-')
plot(toutlin, xoutlin, '*-')
hold off