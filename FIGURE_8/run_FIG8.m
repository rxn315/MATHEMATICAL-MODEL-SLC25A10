clear; clc;

% ---------------- Parameters ----------------
load('MCMC_Result.mat')
Ts_max = posterior_mean(1);  % for succinate
Tm_max = posterior_mean(2);  % for malate

% Michaelis constants (in mM)
Ks_m = posterior_mean(3);  % Succinate
Km_m = posterior_mean(4);    % Malate
Kp_m = posterior_mean(5);    % Phosphate

% Volumes (arbitrary units)
Vm = 1;
Vims = 1/10;

lambda21=1.3373; 
lambda31=0.63160; 

% ---------------- Data / ICs ----------------
load('Palmier_DIC_competition_Exp_Data.mat'); % expects M_with_S
assert(exist('M_with_S','var')==1, 'M_with_S not found in MAT file.');

P_m0 = 10; 
S_m0 = 10;
M_m0 = 10;
P_c0 = 0;
M_c0 = 1./M_with_S(1);
S_c0 = 0.5;

y0 = [M_m0, M_c0, S_m0, S_c0, P_m0, P_c0]; % [Mm, Mims, Sm, Sims, Pm, Pims]

% ---------------- Time span ----------------
% Using minutes here to match your axis labels
tspan = [0 0.05];  % minutes

% ---------------- ODE solve ----------------
[Time, y] = ode15s(@(t,y) slc25a10_model(t, y, Vims, Vm, Ks_m, Km_m, Kp_m, ...
                                         Ts_max, Tm_max, lambda21, lambda31), ...
                   tspan, y0, odeset('RelTol',1e-8,'AbsTol',1e-12));

% ---------------- Extract states ----------------
Mm   = y(:,1);    % Malate (matrix)
Mims = y(:,2);    % Malate (IMS)
Sm   = y(:,3);    % Succinate (matrix)
Sims = y(:,4);    % Succinate (IMS)
Pm   = y(:,5);    % Phosphate (matrix)
Pims = y(:,6);    % Phosphate (IMS)

% ---------------- Gradients (guard tiny denom) ----------------
epsd = 1e-12;
Mal_grad = Mm ./ max(Mims, epsd);
Suc_grad = Sm ./ max(Sims, epsd);
Pho_grad = Pm ./ max(Pims, epsd);


