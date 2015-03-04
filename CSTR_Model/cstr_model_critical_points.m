% Du reactor linearisation
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
            
%% Critical points for different Q            
N = 200;
Qs = linspace(-1272.0, 920.0, N); % move in the range -22kW to 15 kW
xss1 = zeros(N, 2);
xss2 = zeros(N, 2);
xss3 = zeros(N, 2);

Q = Qs(1);
% x(1) = CA
% x(2) = TR
sys = @(t, x) [F/V*CA0 - (F/V + k0*exp(-E/(R*x(2))))*x(1);
                F/V*TA0 - F/V*x(2) - dH/(rho*Cp)*k0*exp(-E/(R*x(2)))*x(1)+Q/(rho*Cp*V)];

guess1 = [0.999, 289.];
guess2 = [0.418, 405.];
guess3 = [0.0115, 486.7];
[xss1t, ~] = fsolve(@(x) sys(1.0, x), guess1, options);
[xss2t, ~] = fsolve(@(x) sys(1.0, x), guess2, options);
[xss3t, ~] = fsolve(@(x) sys(1.0, x), guess3, options);
xss1(1, :) = xss1t;
xss2(1, :) = xss2t;
xss3(1, :) = xss3t;

for k=2:1:N
    Q = Qs(k);
    sys = @(t, x) [F/V*CA0 - (F/V + k0*exp(-E/(R*x(2))))*x(1);
                    F/V*TA0 - F/V*x(2) - dH/(rho*Cp)*k0*exp(-E/(R*x(2)))*x(1)+Q/(rho*Cp*V)];
    guess1 = xss1(k-1, :);
    guess2 = xss2(k-1, :);
    guess3 = xss3(k-1, :);
    [xss1t, ~] = fsolve(@(x) sys(1.0, x), guess1, options);
    [xss2t, ~] = fsolve(@(x) sys(1.0, x), guess2, options);
    [xss3t, ~] = fsolve(@(x) sys(1.0, x), guess3, options);
    xss1(k, :) = xss1t;
    xss2(k, :) = xss2t;
    xss3(k, :) = xss3t;
end

Q = 0.0;
sys = @(t, x) [F/V*CA0 - (F/V + k0*exp(-E/(R*x(2))))*x(1);
                F/V*TA0 - F/V*x(2) - dH/(rho*Cp)*k0*exp(-E/(R*x(2)))*x(1)+Q/(rho*Cp*V)];
guess1 = xss1(round(N/2), :);
guess2 = xss2(round(N/2), :);
guess3 = xss3(round(N/2), :);
[xss1_Q0, ~] = fsolve(@(x) sys(1.0, x), guess1, options);
[xss2_Q0, ~] = fsolve(@(x) sys(1.0, x), guess2, options);
[xss3_Q0, ~] = fsolve(@(x) sys(1.0, x), guess3, options);

%% Plot critical points for different Q
figure(1)
plot(xss1(:,1), xss1(:,2),'b','LineWidth',6)
hold on
plot(xss2(:,1), xss2(:,2),'g','LineWidth',6)
plot(xss3(:,1), xss3(:,2),'r','LineWidth',6)
plot(xss1_Q0(1), xss1_Q0(2),'kx','LineWidth',6,'MarkerSize', 20)
plot(xss2_Q0(1), xss2_Q0(2),'kx','LineWidth',6,'MarkerSize', 20)
plot(xss3_Q0(1), xss3_Q0(2),'kx','LineWidth',6,'MarkerSize', 20)
hold off
xlabel('Concentration [kmol.m^{-3}]','fontsize', 18)
ylabel('Temperature [K]','fontsize', 18)
set(gca,'fontsize', 18);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'cstr_model_stat_points.pdf');

