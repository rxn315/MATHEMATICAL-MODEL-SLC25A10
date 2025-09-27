function createfigure_concentration(Time, Mm, Sm, Pm, Mims, Sims, Pims)
%CREATEFIGURE Plots metabolite concentrations over time
%  Time : vector of time points
%  Mm   : Malate in mitochondria
%  Sm   : Succinate in mitochondria
%  Pm   : Phosphate in mitochondria
%  Mims : Malate in IMS
%  Sims : Succinate in IMS
%  Pims : Phosphate in IMS
%
% Example call:
%   createfigure(Time, Mm, Sm, Pm, Mims, Sims, Pims);

% Create figure
figure('Color',[1 1 1]);

% ---- Row 1 ----
subplot(2,3,1);
plot(Time, Mm, 'k-', 'LineWidth', 2);
title('[M_m]');
xlabel('Time (min)'); ylabel('Concentration (mM)');
set(gca,'FontSize',30,'FontWeight','bold'); grid on; 

subplot(2,3,2);
plot(Time, Sm, 'k-', 'LineWidth', 2);
title('[S_m]');
xlabel('Time (min)'); ylabel('Concentration (mM)');
set(gca,'FontSize',30,'FontWeight','bold'); grid on; 

subplot(2,3,3);
plot(Time, Pm, 'k-', 'LineWidth', 2);
title('[P_m]');
xlabel('Time (min)'); ylabel('Concentration (mM)');
set(gca,'FontSize',30,'FontWeight','bold'); grid on; 

% ---- Row 2 ----
subplot(2,3,4);
plot(Time, Mims, 'k-', 'LineWidth', 2);
title('[M_{ims}]');
xlabel('Time (min)'); ylabel('Concentration (mM)');
set(gca,'FontSize',30,'FontWeight','bold'); grid on; 
subplot(2,3,5);
plot(Time, Sims, 'k-', 'LineWidth', 2);
title('[S_{ims}]');
xlabel('Time (min)'); ylabel('Concentration (mM)');
set(gca,'FontSize',30,'FontWeight','bold'); grid on; 

subplot(2,3,6);
plot(Time, Pims, 'k-', 'LineWidth', 2);
title('[P_{ims}]');
xlabel('Time (min)'); ylabel('Concentration (mM)');
set(gca,'FontSize',30,'FontWeight','bold'); grid on; 

end
