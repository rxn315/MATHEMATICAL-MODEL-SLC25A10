
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

% ------------------------------
% Maximum transport rates
% ------------------------------


R = 8.31e-3; % Gas constant, kJ./(mol·K)
T = 298.15;  % Temperature, K
RT = R .* T;  % RT product, kJ./mol

% Initial Conditions
 P_m0 = 10;
 S_m0 = 10;
 M_m0 = 0;
 P_c0 = 0.5;
 M_c0 = 0;
 S_c0 = 0.1;

 Pm=P_m0;
 Sm=S_m0;
 Mm=M_m0;
 Pims=P_c0;
 Mims=M_c0;
 Sims=S_c0;


% Time span for simulation
tspan = [0 0.06];

% Solve the ODE system
[Time, y] = ode15s(@ (t,y)slc25a10_model(t, y, Vims, Vm, Ks_m, Km_m, Kp_m, Ts_max, Tm_max,lambda21,lambda31), tspan, [M_m0, M_c0, S_m0, S_c0, P_m0, P_c0]);

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
    J_succ = (1/Vm) .* (Ts_max .* (phi2.*Sims - Sm) ./ (Ks_m.*delta1 + Ks_m.*phi2 .* delta2)) ...
           - (1/Vims) .* (Ts_max .* (Sm - phi1.*Sims) ./ (Ks_m.*delta1 + Ks_m.*phi1 .* delta2));

    % J_mal
    J_mal = (1/Vm) .* (Tm_max .* (phi2.*Mims - Mm) ./ (Km_m.*delta1 + Km_m.*phi2 .* delta2)) ...
           - (1/Vims) .* (Tm_max .* (Mm - phi1.*Mims) ./ ( Km_m.*delta1 + Km_m.*phi1 .* delta2));

    % J_pho
    % J_pho = (1/Vims) .* (Tp_max .* (Pm - phi2.*Pims) ./ (Kp_m.*delta1 +  Kp_m.*phi2 .* delta2)) ...
    %        - (1/Vm) .* (Tp_max .* (phi1.*Pims - Pm) ./ ( Kp_m.*delta1 +  Kp_m.*phi1 .* delta2));

    J_S=J_succ;
    J_M=J_mal;
    J_P=-J_succ-J_mal;
    % J_P=J_pho;

   % Overall transporter flux
    J_tot = abs(J_M) + abs(J_S) + abs(J_P);





% Analytical equilibrium concentrations
S_total = S_m0 + S_c0; % Total succinate (mM)
P_total = P_m0 + P_c0; % Total phosphate (mM)
M_total = S_m0 + P_m0; % Total matrix species (mM)
I_total = S_total + P_total - M_total; % Total IMS species (mM)

Succ_m_eq = (S_total .* M_total) ./ (S_total + P_total); % Equilibrium Succ_m
Succ_ims_eq = S_total - Succ_m_eq; % Equilibrium Succ_ims
P_m_eq = M_total - Succ_m_eq; % Equilibrium P_m
P_ims_eq = P_total - M_total + Succ_m_eq; % Equilibrium P_ims

% Calculate Gibbs free energy over time
delta_G = zeros(length(Sm), 1);
for i = 1:length(Sm)
    Succ_mito = Sm(i);
    Mal_mito = Mm(i);
    P_i_mito = Pm(i);
    Succ_cyto = Sims(i);
    Mal_cyto = Mims(i);
    P_i_cyto = Pims(i);
    
    % Avoid log(0) issues
    if P_i_cyto == 0
        P_i_cyto = 1e-6;
    end
    if Succ_mito == 0
        Succ_mito = 1e-6;
    end
    if Succ_cyto == 0
        Succ_cyto = 1e-6;
    end
    if P_i_mito == 0
        P_i_mito = 1e-6;
    end

    delta_G_Succ_mito = -690.44;
    delta_G_Succ_cyto = -690.44;
    delta_G_Mal_mito = -526.16;
    delta_G_Mal_cyto = -526.16;
    delta_G_P_i_mito = -1069.1;
    delta_G_P_i_cyto = -1069.1;

    S = [0, 0, 0, 0];
    s = [-1, -1, 1, 1];
    C = [Succ_cyto, P_i_mito, Succ_mito, P_i_cyto];
    delta_G_prime = [delta_G_Succ_cyto, delta_G_P_i_mito, delta_G_Succ_mito, delta_G_P_i_cyto];

    term1 = RT .* sum(S .* log(C));
    term2 = sum(S .* delta_G_prime);
    term3 = RT .* sum(s .* log(C));
    term4 = sum(s .* delta_G_prime);

    delta_G(i) = term1 + term2 + term3 + term4;
end

% Gibbs free energy as a function of Succ_m
x = linspace(1e-6, min(S_total, M_total) - 1e-6, 1000);
deltaG_func = @(x) RT .* log(((S_total - x) .* (M_total - x)) ./ (x .* (P_total - M_total + x)));
G = deltaG_func(x);


% X axis
X1 = Time(:);

% Build the Y matrices (ensure column vectors)
YMatrix1 = [J_S(:),  J_P(:)];                     % a) fluxes (Jsucc, Jpho)
YMatrix2 = [delta_G(:),  zeros(numel(Time),1)];   % b) ΔG and zero line
YMatrix3 = [Sm(:),   Succ_m_eq*ones(numel(Time),1)];   % c) Sm
YMatrix4 = [Sims(:), Succ_ims_eq*ones(numel(Time),1)]; % d) Sims
YMatrix5 = [Pm(:),   P_m_eq*ones(numel(Time),1)];      % e) Pm
YMatrix6 = [Pims(:), P_ims_eq*ones(numel(Time),1)];    % f) Pims

% Render the figure
createfigure(X1, YMatrix1, YMatrix2, YMatrix3, YMatrix4, YMatrix5, YMatrix6);

