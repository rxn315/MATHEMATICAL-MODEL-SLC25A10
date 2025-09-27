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

% ---------------- Time span (minutes) ----------------
tspan = [0 0.05];  % minutes (match axis labels)

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

% ---------------- Helper terms ----------------
delta1 = 1 + Sm./Ks_m   + Mm./Km_m   + Pm./Kp_m;
delta2 = 1 + Sims./Ks_m + Mims./Km_m + Pims./Kp_m;

epsd = 1e-12;

% phi1
phi1_num = Sm + lambda21.*Mm + lambda31.*Pm;
phi1_den = Sims + lambda21.*Mims + lambda31.*Pims;
phi1 = phi1_num ./ max(phi1_den, epsd);

% phi2 (volume-squared factor on phosphate’s lambda)
vol2    = (Vm / Vims)^2;
phi2_num = Sm + lambda21.*Mm + vol2.*lambda31.*Pm;
phi2_den = Sims + lambda21.*Mims + vol2.*lambda31.*Pims;
phi2 = phi2_num ./ max(phi2_den, epsd);

% ---------------- Fluxes ----------------
den_s_phi2 = Ks_m .* delta1 + Ks_m .* phi2 .* delta2;
den_s_phi1 = Ks_m .* delta1 + Ks_m .* phi1 .* delta2;

den_m_phi2 = Km_m .* delta1 + Km_m .* phi2 .* delta2;
den_m_phi1 = Km_m .* delta1 + Km_m .* phi1 .* delta2;

den_p_phi2 = Kp_m .* delta1 + Kp_m .* phi2 .* delta2;
den_p_phi1 = Kp_m .* delta1 + Kp_m .* phi1 .* delta2;

J_succ  = (Ts_max./Vm)  .* (phi2.*Sims - Sm) ./ max(den_s_phi2, epsd) ...
        - (Ts_max./Vims).* (Sm - phi1.*Sims) ./ max(den_s_phi1, epsd);

J_mal   = (Tm_max./Vm)  .* (phi2.*Mims - Mm) ./ max(den_m_phi2, epsd) ...
        - (Tm_max./Vims).* (Mm - phi1.*Mims) ./ max(den_m_phi1, epsd);

% J_pho   = (Tp_max./Vims).* (Pm - phi2.*Pims) ./ max(den_p_phi2, epsd) ...
%         - (Tp_max./Vm)  .* (phi1.*Pims - Pm) ./ max(den_p_phi1, epsd);
J_pho=-(J_succ+J_mal);
% Total flux (as in your other script: abs-sum)
J_tot = abs(J_succ) + abs(J_mal) + abs(J_pho);

% ---------------- Inputs for your patch+flux figure ----------------
X1 = Time(:);                    % left axis time
YMatrix1 = [J_succ, J_mal, J_pho];  % left axis 3 curves

% Shade efflux (y >= 0) and influx (y <= 0) regions
ymax = max([YMatrix1(:); 0]);
ymin = min([YMatrix1(:); 0]);

XData1 = [X1(1) X1(end) X1(end) X1(1)];  % rectangle in time
YData_efflux = [0 0 ymax ymax];          % above zero (Efflux Region)
YData_influx = [ymin ymin 0 0];          % below zero (Influx Region)

% Right axis (total)
X2 = Time(:);
Y1 = J_tot(:);

