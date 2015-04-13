% LQR system
close all;
clear all;
clc;

offset = [0.48934869384882707;412.1302612302342]; 

A(1,1) = 0.9959087087412309;
A(2,1) = 0.4186578521534139;
A(1,2) = -6.030851990315578e-5;
A(2,2) = 1.010063701982629;

B(1,1) = -2.523369033576597e-9;
B(2,1) = 8.410308376496169e-5;

C = eye(2);
D = 0;
H = [1 0];

% ssmat = [eye(2)-A -B;H*C 0];
% ssvec = [0;0;0]; % want to drive it to zero!
% ssans = ssmat\ssvec

sys = ss(A,B,C,D, 0.1);

Q = zeros(2);
Q(1,1) = 1.0;
R = zeros(1);
R(1,1) = 0.1;

K = dlqr(A,B,Q,R);

% Lets close the loop!
x0 = [0.5;400];
h = 0.1;
tend = 150;
ts = 0:h:tend;
N = length(ts);

xs = zeros(2, N);
xs(:,1) = x0 - offset;
us = zeros(1, N-1);

for t=2:N
    us(:, t-1) = -K*xs(:, t-1); % control law
    xs(:, t) = A*xs(:, t-1) + B*us(:, t-1);
end

xs = xs + repmat(offset, 1, N);

figure(1)
plot(ts, xs(1, :))
figure(2)
plot(ts, xs(2, :))
figure(3)
plot(ts(1:end-1), us)














