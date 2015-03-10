% Runge Kutta for a linearised system
clear all
close all
clc

syms f0 g0 B11 B21 J11 J12 J21 J22 D11 D21 x1 x2 h

F0 = [f0;
      g0];
J = [J11, J12;
     J21, J22];
D = [D11; D21];
B = [B11;B21];

F = @(x) J*x + F0 - D + B;

k1 = F([x1;x2]);
k2 = F([x1;x2] + 1/2*h*k1);
k3 = F([x1;x2] + 1/2*h*k2);
k4 = F([x1;x2] + h*k3);

rk = h/6*(k1+2*k2+2*k3+k4)