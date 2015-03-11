% Runge Kutta for a linearised system
clear all
close all
clc

% syms f0 g0 B11 B21 J11 J12 J21 J22 D11 D21 x1 x2 h
% 
% F0 = [f0;
%       g0];
% J = [J11, J12;
%      J21, J22];
% D = [D11; D21];
% B = [B11;B21];
% 
% F = @(x) J*x + F0 - D + B;
% 
% k1 = F([x1;x2]);
% k2 = F([x1;x2] + 1/2*h*k1);
% k3 = F([x1;x2] + 1/2*h*k2);
% k4 = F([x1;x2] + h*k3);
% 
% rk = simplify((h/6)*(k1+2*k2+2*k3+k4));
% rktop = rk(1,1);
% rkbot = rk(2,1);
% 
% rkcollectedtop = collect(rktop, [x1,x2, B11, B21]);
% rkcollectedbot = collect(rkbot, [x1,x2, B11, B21]);

syms F0 B A D x u h

F = @(x) A*x + B*u +F0 - D;

k1 = F(x);
k2 = F(x + (h/2)*k1);
k3 = F(x + (h/2)*k2);
k4 = F(x + h*k3);

rk = (h/6)*(k1+2*k2+2*k3+k4);
collect(rk, [x, u]) 


