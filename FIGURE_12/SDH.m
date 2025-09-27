 function J_SDH=SDH(Sm,Mm)
% MATLAB script to simulate the flux of succinate dehydrogenase (SDH) | Complex II
% Updated with corrected J_SDH and ksm_eff equations

% Define constants and parameters
R = 8.314; % Gas constant (J/mol·K)
T = 298.15;   % Temperature (K)
DeltaG = -2410; % Gibbs free energy (J/mol)
Keq = exp(-DeltaG / (R * T)); % Equilibrium constant
Vmax = 100; % V_SDH_max (µmol/min/g protein, assumed)
ks_SDH = 1.80; % k^s_SDH (mM)
kq = 140e-3;    % k^q (mM)
kf = 1.80;    % k^f (mM)
kqh = 2.45e-3;   % k^qh (mM)
ko = 0.0315e-3;    % k^o (mM)
km_SDH = 11.7e-3; % k^m_SDH (mM)

% Define concentration ranges
% Sm = linspace(0, 10, 100); % Succinate concentration (0 to 1 mM)
Om = 0;  % Oxaloacetate concentrations (0, 0.01, 0.1 mM)
Qm = 3.75*10^-1;   % Ubiquinone concentration (M, fixed)
Fm = 0.1;   % Fumarate concentration (M, fixed)
QH_m = 1.12; % Ubiquinol concentration (M, fixed)
% Mm = [10, 11, 15];   % Malate concentration (M, fixed)

% Initialize flux array
J_SDH = zeros(length(Sm), length(Om));

% Calculate flux for each [S_m] and [O_m]

    % Compute effective succinate binding constant
    % Corrected model: Oxaloacetate inhibition and Malate stimulation
    ksm_eff = ks_SDH .* ((ko + Om) ./ ko) .* (km_SDH ./ (km_SDH + Mm));
    
    
        % Numerator: Forward - Reverse rates
        forward = (Sm .* Qm) ./ ksm_eff;
        reverse = (Fm .* QH_m) ./ (kf .* kqh .* Keq);
        numerator = Vmax .* (forward - reverse);
        
        % Denominator: Grouped terms as per corrected equation
        denom_part1 = 1 + Sm ./ ksm_eff;
        denom_part2 = (1 + Fm ./ kf) .* (1 + Qm ./ kq + QH_m ./ kqh);
        denominator = denom_part1 + denom_part2;
        
        % Flux calculation
        J_SDH = numerator ./ denominator;
    

 end
% % Plot results
% figure;
% plot(Sm, J_SDH(:,1), 'b-', 'LineWidth', 2, 'DisplayName', '[M_m] = 0.001 mM');
% hold on;
% plot(Sm, J_SDH(:,2), 'r--', 'LineWidth', 2, 'DisplayName', '[M_m] = 0.01 mM');
% plot(Sm, J_SDH(:,3), 'g-.', 'LineWidth', 2, 'DisplayName', '[M_m] = 0.1 mM');
% hold off;
% 
% % Customize plot
% xlabel('Succinate Concentration [S_m] (mM)');
% ylabel('Flux J_{SDH} (\mumol/min/mg protein)');
% %title('SDH Flux vs. Succinate Concentration');
% legend('show');
% grid on;
% set(gca, 'FontSize', 12);
% 
% 
% % % Sm=[0.0108,0.0108,0.0108,0.0108,0.0108,0.0108,0.0108,0.0108,0.0108,0.0108,0.0109,0.011,0.0115,0.0126,0.0149,0.0186,0.024,0.0314,0.0412,0.0529,0.0661,0.0803,0.095,0.1099,0.14,0.1676,0.1922,0.2137,0.2325,0.2628,0.2863,0.3046,0.3191,0.3305,0.3469,0.357,0.3631,0.3668,0.3701,0.3715,0.3722,0.3724];
% % % Time=[0	4.58000000000000e-06	4.58000000000000e-05	0.000417000000000000	0.00380000000000000	0.0338000000000000	0.304000000000000	0.574200000000000	1.15560000000000	1.73710000000000	2.31850000000000	2.89990000000000	3.48140000000000	4.06280000000000	4.64430000000000	5.22570000000000	5.80710000000000	6.38860000000000	7.01260000000000	7.64090000000000	8.27350000000000	8.90610000000000	9.53870000000000	10.1758000000000	11.4967000000000	12.8177000000000	14.1386000000000	15.4596000000000	16.7805000000000	19.3811000000000	21.9817000000000	24.5823000000000	27.1829000000000	29.7835000000000	35.1848000000000	40.5861000000000	45.9873000000000	51.3886000000000	61.1109000000000	70.8332000000000	85.4166000000000	100];
