function J_SCAS=SCAS(Sm,Pm)
% MATLAB script to simulate the flux of SCoA synthetase (SCAS) reaction
% Reaction 6: SCoA_m + GDP_m + Pi_m <=> SUC_m + GTP_m + CoA_m + H_m+

% Define constants
R = 8.314; % Gas constant (J/mol·K)
T = 298.15;   % Temperature (K)
pH_m = 7.4; % Mitochondrial matrix pH
DeltaG0_SCAS = 1.26*10^3; % Gibbs free energy (J/mol, converted from kJ/mol)
Keq_SCAS = exp(-DeltaG0_SCAS / (R * T)) * 10^(pH_m - 7); % Equilibrium constant

% Kinetic parameters
Vmaxf_SCAS = 100; % V_maxf_SCAS (µmol/min/g protein, assumed)
KA = 12.5e-3;  % SCoA binding constant (mM)
KB = 5.70e-3;  % GDP binding constant (mM)
KC = 2.10;  % Pi binding constant (mM)
KD = 451e-3;   % SUC binding constant (mM)
KE = 23.1e-3;  % GTP binding constant (mM)
KF = 17.7e-3;  % CoA binding constant (mM)

% Fixed mitochondrial concentrations
GDP_m = 2.5e-2;   % [GDP_m] (mM)
CoA_m = 10;      % [CoA_m] (mM)
SCoA_m = 10000;   % [SCoA_m] (mM)
GTP_m = 9.75e-1;  % [GTP_m] (mM)

% Variable concentrations
% Sm = linspace(0, 20, 100); % SCoA concentration (0 to 20 mM)
% Pm = [1,5,10];     % CoA concentration (5, 10, 15 mM)

% Initialize flux array
J_SCAS = zeros(length(Sm), length(Pm));

% Calculate flux for each [SCoA_m] and [CoA_m]

        % Numerator: Forward - Reverse rates
        forward = (SCoA_m .* GDP_m .* Pm) ./ (KA .* KB .* KC);
        reverse = (Sm .* GTP_m .* CoA_m) ./ (KD .* KE .* KF .* Keq_SCAS);
        numerator = Vmaxf_SCAS .* (forward - reverse);
        
        % Denominator
        denom_part1 = 1 + SCoA_m./KA + Sm./KD + CoA_m./KF + (Sm .* CoA_m)./(KD .* KF);
        denom_part2 = 1 + GDP_m./KB + Pm./KC + GTP_m./KE + (GDP_m .* Pm)./(KB .* KC);
        denominator = denom_part1 .* denom_part2;
        
        % Flux calculation
        J_SCAS = numerator ./ denominator;
   
end
% % Plot results
% figure;
% plot(Sm, J_SCAS(:,1), 'b-', 'LineWidth', 2, 'DisplayName', '[P_m] = 1 mM');
% hold on;
% plot(Sm, J_SCAS(:,2), 'r--', 'LineWidth', 2, 'DisplayName', '[P_m] = 5 mM');
% plot(Sm, J_SCAS(:,3), 'g-.', 'LineWidth', 2, 'DisplayName', '[P_m] = 10 mM');
% hold off;
% 
% % Customize plot
% xlabel('Succinate Concentration [S_m] (mM)');
% ylabel('Flux J_{SCAS} (\mumol/min/mg protein)');
% title('SCAS Flux vs. Succinate Concentration');
% legend('show');
% grid on;
% set(gca, 'FontSize', 12);