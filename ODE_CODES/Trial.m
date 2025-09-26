clear; clc;

% ---------------- Parameters ----------------
load('MCMC_Result.mat')
Ts_max = posterior_mean(1);  % for succinate
Tm_max = posterior_mean(2);  % for malate
Tp_max =100;

% Michaelis constants (in mM)
Ks_m = posterior_mean(3);  % Succinate
Km_m = posterior_mean(4);    % Malate
Kp_m = posterior_mean(5);    % Phosphate

% Volumes (arbitrary units)
Vm = 1;
Vims = 1/10;

lambda21=1.3373; 
lambda31=0.63160; 




% Experimental data (replace with actual values)
load('Palmier_DIC_competition_Exp_Data.mat')
 P_m0 = 10;
S_m0 = 10;
M_m0 = 10;
P_c0 = 0;
M_c0 = 1./M_with_S(1);
S_c0 = 0.5;



% Time span for simulation
tspan = [0 0.1];

% Solve the ODE system
[Time, y] = ode15s(@ (t,y)slc25a10_model(t, y, Vims, Vm, Ks_m, Km_m, Kp_m, Ts_max, Tm_max, Tp_max,lambda21,lambda31), tspan, [M_m0, M_c0, S_m0, S_c0, P_m0, P_c0]);

% Extract results
Mm = y(:, 1); % Malate in mitochondria
Mims = y(:, 2); % Malate in cytoplasm
Sm = y(:, 3); % Succinate in mitochondria
Sims = y(:, 4); % Succinate in cytoplasm
Pm = y(:, 5); % Phosphate in mitochondria
Pims = y(:, 6); % Phosphate in cytoplasm


% -----------------------------
    % Delta Terms
    % -----------------------------
    delta1 = 1 + Sm./Ks_m + Mm./Km_m + Pm./Kp_m;
    delta2 = 1 + Sims./Ks_m + Mims./Km_m + Pims./Kp_m;

    % -----------------------------
    % Phi Terms
    % -----------------------------
    % phi1
    phi1_num = Sm + lambda21.*Mm + lambda31.*Pm;
    phi1_den = Sims + lambda21.*Mims + lambda31.*Pims;
    phi1 = phi1_num ./ phi1_den;

    % phi2
    vol2 = (Vm / Vims)^2;
    phi2_num = Sm + lambda21.*Mm + vol2.*lambda31.*Pm;
    phi2_den = Sims + lambda21.*Mims + vol2.*lambda31.*Pims;
    phi2 = phi2_num ./ phi2_den;
    % -----------------------------
    % Flux Equations
    % -----------------------------
    % J_succ
    J_succ = (1./Vm) .* (Ts_max .* (phi2.*Sims - Sm) ./ (Ks_m .* delta1 + Ks_m .* phi2 .* delta2)) ...
           - (1./Vims) .* (Ts_max .* (Sm - phi1.*Sims) ./ (Ks_m .* delta1 + Ks_m .* phi1 .* delta2));

    % J_mal
    J_mal = (1./Vm) .* (Tm_max .* (phi2.*Mims - Mm) ./ (Km_m .* delta1 + Km_m .* phi2 .* delta2)) ...
           - (1./Vims) .* (Tm_max .* (Mm - phi1.*Mims) ./ (Km_m .* delta1 + Km_m .* phi1 .* delta2));

    % J_pho
    J_pho = (1./Vims) .* (Tp_max .* (Pm - phi2.*Pims) ./ (Kp_m .* delta1 + Kp_m .* phi2 .* delta2)) ...
           - (1./Vm) .* (Tp_max .* (phi1.*Pims - Pm) ./ (Kp_m .* delta1 + Kp_m .* phi1 .* delta2));


 