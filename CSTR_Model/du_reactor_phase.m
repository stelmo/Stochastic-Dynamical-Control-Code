% Du reactor phase portrait
clear all
close all
clc

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

f = @(t, x) [-phi*x(1)*k(x(2))+q*(x1f-x(1)); 
             beta*phi*x(1)*k(x(2))-(q+delta)*x(2)+delta*u+q*x2f]; % For ODE solving
F = @(x) [-phi*x(1)*k(x(2))+q*(x1f-x(1)); 
          beta*phi*x(1)*k(x(2))-(q+delta)*x(2)+delta*u+q*x2f]; % For SS solving

%% Solve for the steady states
ss1 = [0.85; 0.88];
ss2 = [0.55; 2.75];
ss3 = [0.23; 4.71];
[xss1, fval] = fsolve(F,ss1);
[xss2, fval] = fsolve(F,ss2);
[xss3, fval] = fsolve(F,ss3);

%% Create a phase portrait
x1 = linspace(0.0, 1, 10);
x2 = linspace(0.0, 6, 10); 
% [tout, xs] = ode45(f, [0 100], [0.24; 4.8]);

[x1s, x2s] = meshgrid(x1, x2);
u = zeros(size(x1s));
v = zeros(size(x2s));

t=0.0;
for i=1:numel(x1s)
   xprime = f(t, [x1s(i); x2s(i)]);
   u(i) = xprime(1);
   v(i) = xprime(2);
end

%% Plotting
figure(1)
quiver(x1s, x2s, u, v,'r')
% plot(xs(:,1), xs(:,2), 'r')
hold on
plot(xss1(1), xss1(2), 'bx');
% plot(xss2(1), xss2(2), 'bx');
% plot(xss3(1), xss3(2), 'bx');
hold off
xlabel('x_1')
ylabel('x_2')
% axis tight equal;