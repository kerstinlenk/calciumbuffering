%% Author of the code: Kerstin Lenk
%  Seleton of the code produced with Perplexity, April 2026
%  Related publication: 
%  Calcium buffering in astrocytes and its relevance for experimental data interpretation and computational modelling
%  Kerstin Lenk*, Andre Zeug, Franziska E. Müller

%% Calcium buffer schematic figure
% Panel A: same Kd, different kon/koff
% Panel B: different Kd vs EC50 combinations and their effect on Ca2+ transient

clear; close all; clc;

%% ------------------------------------------------------------
%% Panel A: Same Kd, different kinetic speeds

%% Parameters
t_span = [0 15];                 % simulate 20 seconds
K_d = 50;                        % uM, low affinity (same for both buffers)
B_total = 30;                    % µM, total buffer concentration
ca_baseline = 0.084;             % µM, baseline Ca2+; Shigetomi et al. 2016

%% Buffer kinetics (same Kd, different kon)
% Fast buffer
kon_fast  = 50;                  % µM^-1 s^-1
koff_fast = K_d * kon_fast;

% Slow buffer
kon_slow  = 0.1;                 % µM^-1 s^-1
koff_slow = K_d * kon_slow;

%% Oscillatory Ca2+ input
A = 1.0;                         % µM, amplitude
tau_r = 0.03;                    % s, rise time
tau_d = 0.45;                    % s, decay time
f = 0.2;                         % Hz

ca_input = @(t) ...
    max(0, A*(exp(-(mod(t,1/f))/tau_d) - exp(-(mod(t,1/f))/tau_r)));

%% ODE systems
% 1. No buffer
ode_none = @(t,y) ca_input(t) - 2*(y(1)-ca_baseline);

% 2. Fast buffer
ode_fast = @(t,y) [ca_input(t) - 4*(y(1)-ca_baseline) ...
    - (kon_fast*y(1)*(B_total-y(2)) - koff_fast*y(2));
    kon_fast*y(1)*(B_total-y(2)) - koff_fast*y(2)
];

% 3. Slow buffer
ode_slow = @(t,y) [ca_input(t) - 2*(y(1)-ca_baseline) ...
    - (kon_slow*y(1)*(B_total-y(2)) - koff_slow*y(2));
    kon_slow*y(1)*(B_total-y(2)) - koff_slow*y(2)
];

%% Solve
[t1,y1] = ode45(ode_none, t_span, ca_baseline);
[t2,y2] = ode45(ode_fast, t_span, [ca_baseline; 0]);
[t3,y3] = ode45(ode_slow, t_span, [ca_baseline; 0]);

%% ------------------------------------------------------------
%% Panel B: Different Kd vs EC50 combinations
%  Two indicator scenarios:
%  1) strong binder but fluorescence responds later (Kd << EC50)
%  2) weaker binder but fluorescence responds earlier (Kd > EC50)
% 
%  Binding changes the Ca2+ transient; fluorescence is computed with a Hill curve:
%  Fnorm = Ca^n / (EC50^n + Ca^n)

nH     = 2.5;    % Hill coefficient for fluorescence readout

ind(1).name = 'Kd << EC50';
ind(1).Kd   = 0.15;      % uM
ind(1).EC50 = 0.80;      % uM
ind(1).kon  = 35;        % 1/(uM*s)
ind(1).koff = ind(1).Kd * ind(1).kon;

ind(2).name = 'Kd > EC50';
ind(2).Kd   = 0.90;      % uM
ind(2).EC50 = 0.35;      % uM
ind(2).kon  = 20;        % 1/(uM*s)
ind(2).koff = ind(2).Kd * ind(2).kon;

colorsB = [0.55 0.10 0.70;
           0.00 0.60 0.60];

%% Dose-response curves for panel B
Ca_scan = logspace(-2, 1, 500); % 0.01 to 10 uM
theta1 = Ca_scan ./ (ind(1).Kd + Ca_scan); % occupancy proxy for Kd
theta2 = Ca_scan ./ (ind(2).Kd + Ca_scan);

Fscan1 = (Ca_scan.^nH) ./ (ind(1).EC50^nH + Ca_scan.^nH);
Fscan2 = (Ca_scan.^nH) ./ (ind(2).EC50^nH + Ca_scan.^nH);

%% ------------------------------------------------------------
%% Plot
figure('Color','w','Position',[100 100 1500 600]);

tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% A: free Ca transient, same Kd different kinetics
nexttile;
plot(t1, y1, 'k--', 'LineWidth', 1.5); hold on;
plot(t2, y2(:,1), 'r', 'LineWidth', 2.5);
plot(t3, y3(:,1), 'b', 'LineWidth', 2.5);
xlabel('Time (s)', 'FontSize', 16);
ylabel('Free Ca^{2+} (\muM)', 'FontSize', 16);
title('(A) Same K_d, different k_{on}/k_{off}', 'FontSize', 16);
legend('Native transient', 'Fast buffer', 'Slow buffer', 'Location','southeast', 'FontSize', 12);
box off;

% B: Kd occupancy vs EC50 fluorescence curves
nexttile;
% yyaxis left
semilogx(Ca_scan, theta1, '-', 'Color', colorsB(1,:), 'LineWidth', 2.5); hold on;
xline(ind(1).Kd,   ':', 'Color', colorsB(1,:), 'LineWidth', 1.);
xline(ind(1).EC50, '--', 'Color', colorsB(1,:), 'LineWidth', 1.);

semilogx(Ca_scan, theta2, '-', 'Color', colorsB(2,:), 'LineWidth', 2.5);
ylabel('Binding occupancy', 'FontSize', 16);
ylim([0 1]);
xline(ind(2).Kd,   ':', 'Color', colorsB(2,:), 'LineWidth', 1.);
xline(ind(2).EC50, '--', 'Color', colorsB(2,:), 'LineWidth', 1.);
hold off;

xlabel('[Ca^{2+}] (\muM)', 'FontSize', 16);
title('(B) K_d versus EC_{50}', 'FontSize', 16);
legend({'K_d<<EC_{50} case: Binding occ.', ...
    'K_d<<EC_{50} case: K_d','K_d<<EC_{50} case: EC_{50}', ...
    'K_d>EC_{50} case: Binding occ.',...
    'K_d>EC_{50} case: K_d', 'K_d>EC_{50} case: EC_{50}'}, ...
       'Location','southeast', 'FontSize', 12);
box off; 

%% Optional export
exportgraphics(gcf, 'buffer_Kd_EC50_schematic.png', 'Resolution', 600);