% ===========================================
% Local Sensitivity Analysis for J_mal (reparameterized)
% Params: [Ts_max, Tm_max, Ks_m, Km_m, Kp_m, lambda21, lambda31]
% ===========================================

clear; clc;

% Load experimental competition data
load('Palmier_DIC_Competition_Exp_Data.mat');  % provides M_with_S, etc.

% ---------- Fixed volumes ----------
Vm   = 1;
Vims = 1/10;

% ---------- Fixed phosphate transport capacity (not estimated) ----------
TP_FIXED = 1;   % <-- set to your adopted literature value and units

% ---------- Parameter names (7 parameters) ----------
param_names = { ...
    'T^{s}_{max}', 'T^{m}_{max}', ...
    'K^{s}_{m}',   'K^{m}_{m}',   'K^{p}_{m}', ...
    '\lambda_{21}', '\lambda_{31}'};
n_params = numel(param_names);

% ---------- Baseline parameters (edit to your calibrated values) ----------
% Order: [Ts_max, Tm_max, Ks_m, Km_m, Kp_m, lambda21, lambda31]
params0 = [70, 70, 1.17, 0.23, 0.93, 1.0, 1.0];

% ---------- Finite-difference perturbation ----------
perturbation = 0.01;  % ±1%

% ---------- Preallocate for 6 experimental conditions ----------
LSC_MATRIX = zeros(6, n_params);  % LSCs for J_mal per condition

% ---------- Loop over 6 experimental conditions ----------
for j = 1
    % Initial conditions for condition j (match your dataset)
    P_m0 = 10;
    S_m0 = 10;
    M_m0 = 10;
    P_c0 = 0;
    M_c0 = 1 / M_with_S(j);
    S_c0 = 0.5;
    y0 = [M_m0, M_c0, S_m0, S_c0, P_m0, P_c0]; % [Mm, Mc, Sm, Sc, Pm, Pc]

    % Compute baseline malate flux
    [~, J_mal_base, ~, ~] = compute_flux(y0, params0, Vm, Vims, TP_FIXED);

    % Local Sensitivity Coefficients for J_mal
    LSC_vec = zeros(1, n_params);

    for i = 1:n_params
        % Perturb parameter i by ±1%
        p_plus  = params0; p_plus(i)  = params0(i) * (1 + perturbation);
        p_minus = params0; p_minus(i) = params0(i) * (1 - perturbation);

        % Perturbed fluxes
        [~, J_mal_plus,  ~, ~] = compute_flux(y0, p_plus,  Vm, Vims, TP_FIXED);
        [~, J_mal_minus, ~, ~] = compute_flux(y0, p_minus, Vm, Vims, TP_FIXED);

        % Central finite difference derivative dJ/dp_i
        dJdpi = (J_mal_plus - J_mal_minus) / (2 * perturbation * params0(i));

        % Normalized local sensitivity coefficient for J_mal
        LSC_vec(i) = (params0(i) / J_mal_base) * dJdpi;
    end

    % Store this condition’s LSCs
    LSC_MATRIX(j, :) = LSC_vec;
end

% ---------- Mean LSCs across all conditions ----------
mean_LSCs_mal = mean(LSC_MATRIX, 1);
save('LSCs_mal.mat','mean_LSCs_mal');

% ---------- Plot horizontal bar chart ----------
figure;
barh(mean_LSCs_mal, 'FaceColor', [0.2, 0.4, 0.6]);
set(gca, 'YTick', 1:n_params, 'YTickLabel', param_names, 'FontSize', 14);
xlabel('Normalized Sensitivity (LSC) for J_{mal}', 'FontSize', 16);
ylabel('Parameters', 'FontSize', 16);
title('Mean Local Sensitivity Coefficients for Malate Flux', 'FontSize', 18);
grid on;

% ===========================
% Flux computation function
% ===========================
function [J_succ, J_mal, J_pho, J_total] = compute_flux(y, params, Vm, Vims, TP_FIXED)
    % States
    Mm  = y(1); Mims = y(2); Sm  = y(3); Sims = y(4); Pm = y(5); Pims = y(6);

    % Parameters (reparameterized set)
    Ts_max   = params(1);
    Tm_max   = params(2);
    Ks_m     = params(3);   % succinate
    Km_m     = params(4);   % malate
    Kp_m     = params(5);   % phosphate
    lambda21 = params(6);   % = lambda2 / lambda1
    lambda31 = params(7);   % = lambda3 / lambda1

    % ----- Delta terms -----
    delta1 = 1 + Sm/Ks_m + Mm/Km_m + Pm/Kp_m;
    delta2 = 1 + Sims/Ks_m + Mims/Km_m + Pims/Kp_m;

    % ----- Phi terms (using lambda ratios) -----
    phi1 = (Sm + lambda21*Mm + lambda31*Pm) / ...
           (Sims + lambda21*Mims + lambda31*Pims);

    vol2 = (Vm / Vims)^2;  % scaling for phosphate in phi2
    phi2 = (Sm + lambda21*Mm + vol2*lambda31*Pm) / ...
           (Sims + lambda21*Mims + vol2*lambda31*Pims);

    % ----- Fluxes (full forward-minus-reverse forms) -----
    % Succinate
    J_succ = (1/Vm)  * (Ts_max * (phi2*Sims - Sm) / (Ks_m * (delta1 + phi2 * delta2))) ...
           - (1/Vims)* (Ts_max * (Sm - phi1*Sims) / (Ks_m * (delta1 + phi1 * delta2)));

    % Malate
    J_mal  = (1/Vm)  * (Tm_max * (phi2*Mims - Mm) / (Km_m * (delta1 + phi2 * delta2))) ...
           - (1/Vims)* (Tm_max * (Mm - phi1*Mims) / (Km_m * (delta1 + phi1 * delta2)));

    % Phosphate (TP_FIXED is constant, not estimated)
    J_pho  = (1/Vims)* (TP_FIXED * (Pm - phi2*Pims) / (Kp_m * (delta1 + phi2 * delta2))) ...
           - (1/Vm)  * (TP_FIXED * (phi1*Pims - Pm) / (Kp_m * (delta1 + phi1 * delta2)));

    % Total magnitude (optional)
    J_total = abs(J_succ) + abs(J_mal) + abs(J_pho);
end
