%% Initial Fluxes under Two Conditions + Shaded Regions
% Condition A: Vm varies (1:0.1:10), Vims = Vm/10       → black line, right shaded band (10^0..10^1)
% Condition B: Vm = 1 fixed, Vims varies (0.1..10)      → red dashed, left shaded band (10^-1..10^0)

clear; clc;

% ---------------- Parameters ----------------
load('MCMC_Result.mat')
Ts_max = posterior_mean(1);  % for succinate
Tm_max = posterior_mean(2);  % for malate

% Michaelis constants (in mM)
Ks_m = posterior_mean(3);  % Succinate
Km_m = posterior_mean(4);    % Malate
Kp_m = posterior_mean(5);    % Phosphate

lambda21=1.3373; 
lambda31=0.63160; 

% ---------- Initial concentrations ----------
load('Palmier_DIC_competition_Exp_Data.mat'); % for M_with_S
P_m0 = 10;
S_m0 = 10;
M_m0 = 10;
P_c0 = 0.1;
M_c0 = 1./M_with_S(1);
S_c0 = 0.5;
y0 = [M_m0, M_c0, S_m0, S_c0, P_m0, P_c0]; % [Mm, Mims, Sm, Sims, Pm, Pims]

% ---------- Condition A: Vm varies, Vims = Vm/10 ----------
Vm_A   = 1:0.1:10;
Vims_A = Vm_A/10;
nA = numel(Vm_A);

J_succ_A = zeros(nA,1); J_mal_A = zeros(nA,1); J_pho_A = zeros(nA,1); J_tot_A = zeros(nA,1);

for k = 1:nA
    [Js, Jm, Jp] = fluxes_at_state(y0, Vm_A(k), Vims_A(k), ...
        Ts_max,Tm_max,Ks_m,Km_m,Kp_m,lambda21,lambda31);
    J_succ_A(k) = Js; 
    J_mal_A(k)  = Jm; 
    J_pho_A(k)  = Jp; 
    J_tot_A(k)  = abs(Js)+abs(Jm)+abs(Jp);
end

% ---------- Condition B: Vm fixed = 1, Vims varies ----------
% Vm_B   = ones(1,100);
% Vims_B = logspace(-1, 1, 100);   % 0.1 .. 10 on log scale
% nB = numel(Vims_B);


Vm_B   = ones(1,100);
Vims_B = linspace(0.1,1,100);
nB = numel(Vims_B);

J_succ_B = zeros(nB,1); J_mal_B = zeros(nB,1); J_pho_B = zeros(nB,1); J_tot_B = zeros(nB,1);

for k = 1:nB
    [Js, Jm, Jp] = fluxes_at_state(y0, Vm_B(k), Vims_B(k), ...
        Ts_max,Tm_max,Ks_m,Km_m,Kp_m,lambda21,lambda31);
    J_succ_B(k) = Js; 
    J_mal_B(k)  = Jm; 
    J_pho_B(k)  = Jp; 
    J_tot_B(k)  = abs(Js)+abs(Jm)+abs(Jp);
end

%% ---------- Plot results with shaded regions ----------
figure('Color',[1 1 1]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

xL = [1e-1 1e1];      % display 10^-1 to 10^1
alphaVal = 0.10;      % transparency for shading

make_panel(Vm_A,  J_succ_A, Vims_B, J_succ_B, 'J_{succ}(t=0)', 'Initial Succinate Flux', xL, alphaVal);
make_panel(Vm_A,  J_mal_A,  Vims_B, J_mal_B,  'J_{mal}(t=0)',  'Initial Malate Flux',   xL, alphaVal);
make_panel(Vm_A,  J_pho_A,  Vims_B, J_pho_B,  'J_{pho}(t=0)',  'Initial Phosphate Flux',xL, alphaVal);
make_panel(Vm_A,  J_tot_A,  Vims_B, J_tot_B,  'J_{tot}(t=0)',  'Initial Total Transport Activity', xL, alphaVal);

% One legend for all panels (grab handles from last axes)
ax = gca;
legend(ax, {'Vm = 1, Vims↑','Vm↑, Vims = Vm/10'}, ...
    'Location','southoutside','Orientation','horizontal');

%% ===== Helper: flux calculation =====
function [J_succ, J_mal, J_pho] = fluxes_at_state(y, Vm, Vims, ...
        Ts_max,Tm_max,Ks_m,Km_m,Kp_m,lambda21,lambda31)

    Mm   = y(1); Mims = y(2);
    Sm   = y(3); Sims = y(4);
    Pm   = y(5); Pims = y(6);

    % Delta terms
    delta1 = 1 + Sm./Ks_m + Mm./Km_m + Pm./Kp_m;
    delta2 = 1 + Sims./Ks_m + Mims./Km_m + Pims./Kp_m;

    % Phi terms
    phi1_num = Sm + lambda21.*Mm + lambda31.*Pm;
    phi1_den = Sims + lambda21.*Mims + lambda31.*Pims;
    phi1     = phi1_num ./ max(phi1_den, eps);

    vol2     = (Vm / Vims)^2;
    phi2_num = Sm + lambda21.*Mm + vol2.*lambda31.*Pm;
    phi2_den = Sims + lambda21.*Mims + vol2.*lambda31.*Pims;
    phi2     = phi2_num ./ max(phi2_den, eps);

    % Fluxes
    J_succ = (1./Vm)  .* (Ts_max .* (phi2.*Sims - Sm) ./ (Ks_m .* (delta1 + phi2 .* delta2))) ...
           - (1./Vims).* (Ts_max .* (Sm - phi1.*Sims) ./ (Ks_m .* (delta1 + phi1 .* delta2)));

    J_mal  = (1./Vm)  .* (Tm_max .* (phi2.*Mims - Mm) ./ (Km_m .* (delta1 + phi2 .* delta2))) ...
           - (1./Vims).* (Tm_max .* (Mm - phi1.*Mims) ./ (Km_m .* (delta1 + phi1 .* delta2)));

    % J_pho  = (1./Vims).* (Tp_max .* (Pm - phi2.*Pims) ./ (Kp_m .* (delta1 + phi2 .* delta2))) ...
    %        - (1./Vm)  .* (Tp_max .* (phi1.*Pims - Pm) ./ (Kp_m .* (delta1 + phi1 .* delta2)));

    % Electrically neutral alternative:
     J_pho = -(J_succ+J_mal);
end

%% ===== Helper: single panel with shading =====
function make_panel(xA,yA,xB,yB,ylab,ttl,xL,alphaVal)
    nexttile; ax = gca; hold(ax,'on');  grid(ax,'on');
    set(ax,'XScale','log'); xlim(ax,xL);

    % 1) Plot curves first (autoscale y)
    plot(ax, xB, yB, 'r--', 'LineWidth', 2, 'DisplayName','Vm = 1, Vims↑');    % left region
    plot(ax, xA, yA, 'k-',  'LineWidth', 2, 'DisplayName','Vm↑, Vims = Vm/10'); % right region

    % 2) Read y-limits and shade regions
    yl = ylim(ax);

    % Left shaded band: 10^-1 .. 10^0 for Vims sweep
    xb = [1e-1 1e0 1e0 1e-1];
    yb = [yl(1) yl(1) yl(2) yl(2)];
    pB = fill(ax, xb, yb, [1 0 0], 'FaceAlpha', alphaVal, 'EdgeColor','none', 'HandleVisibility','off');

    % Right shaded band: 10^0 .. 10^1 for Vm sweep
    xa = [1e0 1e1 1e1 1e0];
    ya = [yl(1) yl(1) yl(2) yl(2)];
    pA = fill(ax, xa, ya, [0 0 0], 'FaceAlpha', alphaVal, 'EdgeColor','none', 'HandleVisibility','off');

    % 3) Send patches behind lines
    uistack(pB, 'bottom');
    uistack(pA, 'bottom');

    xlabel(ax,'Volume (Vm or Vims)'); ylabel(ax,ylab); title(ax,ttl);
end
