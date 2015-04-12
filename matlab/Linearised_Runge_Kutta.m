% Runge Kutta for a linearised system
clear all
close all
clc

syms B A x u h

F = @(x) A*x + B*u;

k1 = F(x);
k2 = F(x + (h/2)*k1);
k3 = F(x + (h/2)*k2);
k4 = F(x + h*k3);

rk = (h/6)*(k1+2*k2+2*k3+k4);
collect(rk, [x, u]) 


