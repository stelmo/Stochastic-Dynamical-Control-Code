% Wood Berry column
clear all
close all
clc
format long
s = tf('s');


Gc = [12.8/(16.7*s+1) -18.9/(21.0*s+1); 6.6/(10.9*s+1) -19.4/(14.4*s+1)];

[Ac, Bc, Cc, Dc] = ssdata(Gc); % continuous time

Gd = c2d(Gc, 0.5);

[Ad, Bd, Cd, Dd] = ssdata(Gd)

dlmwrite('WB_A.csv', Ad, 'precision', '%.10f')
dlmwrite('WB_B.csv', Bd, 'precision', '%.10f')
dlmwrite('WB_C.csv', Cd, 'precision', '%.10f')
dlmwrite('WB_D.csv', Dd, 'precision', '%.10f')