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





% Experimental data (replace with actual values)
load('Palmier_DIC_competition_Exp_Data.mat')
 P_m0 = 10;
S_m0 = 10;
M_m0 = 10;
P_c0 = 0.1;
M_c0 = 1./M_with_S(1);
S_c0 = 0.5;



% Time span for simulation
tspan = [0 1];

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
    J_succ = (1./Vm) .* (Ts_max .* (phi2.*Sims - Sm) ./ (Ks_m .* delta1 + Ks_m .* phi2 .* delta2)) ...
           - (1./Vims) .* (Ts_max .* (Sm - phi1.*Sims) ./ (Ks_m .* delta1 + Ks_m .* phi1 .* delta2));

    % J_mal
    J_mal = (1./Vm) .* (Tm_max .* (phi2.*Mims - Mm) ./ (Km_m .* delta1 + Km_m .* phi2 .* delta2)) ...
           - (1./Vims) .* (Tm_max .* (Mm - phi1.*Mims) ./ (Km_m .* delta1 + Km_m .* phi1 .* delta2));

    % J_pho
    % J_pho = (1./Vims) .* (Tp_max .* (Pm - phi2.*Pims) ./ (Kp_m .* delta1 + Kp_m .* phi2 .* delta2)) ...
    %        - (1./Vm) .* (Tp_max .* (phi1.*Pims - Pm) ./ (Kp_m .* delta1 + Kp_m .* phi1 .* delta2));


    J_S=J_succ;
    J_M=J_mal;
    % J_P=(1-exp(-Pm)).*(-J_succ-J_mal);
    J_P=(-J_succ-J_mal);
    %J_P=J_pho;

   % Overall transporter flux
    J_tot = abs(J_M) + abs(J_S) + abs(J_P);
     

    

hold on 
% Plot the results
figure;
subplot(2, 3, 1);
plot(Time, Mm);
title('[M_m]');
xlabel('Time (min)');
ylabel('Concentration (mM)');
% legend('M_m (mitochondria)', 'M_c (cytoplasm)');

subplot(2, 3, 4);
plot(Time, Mims);
title('[M_{ims}]');
xlabel('Time (min)');
ylabel('Concentration (mM)');
% legend('M_m (mitochondria)', 'M_c (cytoplasm)');




subplot(2, 3, 2);
plot(Time, Sm);
title('[S_m]');
xlabel('Time (min)');
ylabel('Concentration (mM)');
% legend('S_m (mitochondria)', 'S_c (cytoplasm)');

subplot(2, 3, 5);
plot(Time, Sims);
title('[S_{ims}]');
xlabel('Time (min)');
ylabel('Concentration (mM)');
% legend('S_m (mitochondria)', 'S_c (cytoplasm)');

subplot(2, 3, 3);
plot(Time, Pm);
title('[P_m]');
xlabel('Time (min)');
ylabel('Concentration (mM)');
% legend('P_m (mitochondria)', 'P_c (cytoplasm)');

subplot(2, 3, 6);
plot(Time, Pims);
title('[P_{ims}]');
xlabel('Time (min)');
ylabel('Concentration (mM)');
% legend('P_m (mitochondria)', 'P_c (cytoplasm)');


figure
fill([Time(1), Time(end), Time(end), Time(1)], [-210, -210, 0, 0], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
hold on 
fill([Time(1), Time(end), Time(end), Time(1)], [0, 0, 210, 210], [1 0.6 0], 'FaceAlpha', 0.1, 'EdgeColor', 'none');

plot(Time, J_S);
hold on 
plot(Time, J_M);
hold on
plot(Time, J_P);
xlabel('Time (min)');
ylabel('Flux (\mumol/min/g)');
legend('Efflux Region','Influx Region','J_{succ}','J_{mal}','J_{pho}');

% figure
% plot(Time, SDH(Sm,Mm));
% xlabel('Time (min)');
% ylabel('Flux (\mumol/min/g)');
% hold on
% plot(Time, SCAS(Sm,Pm));
% xlabel('Time (min)');
% ylabel('Flux (\mumol/min/g)');
% legend('J_{SCAS}');

% figure
% subplot(1,3,1)
% plot3(SDH(Sm,Mm),SCAS(Sm,Pm),J_S)
% subplot(1,3,2)
% plot3(SDH(Sm,Mm),SCAS(Sm,Pm),J_M)
% subplot(1,3,3)
% plot3(SDH(Sm,Mm),SCAS(Sm,Pm),J_P)
% 
% figure
% subplot(3,2,1)
% plot(SDH(Sm,Mm),J_S)
% subplot(3,2,2)
% plot(SCAS(Sm,Pm),J_S)
% subplot(3,2,3)
% plot(SDH(Sm,Mm),J_M)
% subplot(3,2,4)
% plot(SCAS(Sm,Pm),J_M)
% subplot(3,2,5)
% plot(SDH(Sm,Mm),J_P)
% subplot(3,2,6)
% plot(SCAS(Sm,Pm),J_P)
% 
% figure
% subplot(3,2,1)
% plot(SDH(Sm,Mm),Sm)
% hold on
% plot(SCAS(Sm,Pm),Sm)
% subplot(3,2,2)
% plot(SDH(Sm,Mm),Sc)
% hold on
% plot(SCAS(Sm,Pm),Sc)
% subplot(3,2,3)
% plot(SDH(Sm,Mm),Mm)
% hold on 
% plot(SCAS(Sm,Pm),Mm)
% subplot(3,2,4)
% plot(SDH(Sm,Mm),Mc)
% hold on 
% plot(SCAS(Sm,Pm),Mc)
% subplot(3,2,5)
% plot(SDH(Sm,Mm),Pm)
% hold on
% plot(SCAS(Sm,Pm),Pm)
% subplot(3,2,6)
% plot(SDH(Sm,Mm),Pc)
% subplot(3,2,6)
% plot(SCAS(Sm,Pm),Pc)
% 
figure
plot(Time,J_tot)
hold on 
% plot(Time,zeros(1,length(SDH(Sm,Mm))))
plot(Time,SDH(Sm,Mm))
hold on 
plot(Time,SCAS(Sm,Pm))
xlabel('Time (min)');
ylabel('Flux (\mumol/min/g)');
legend('J_{tot}','J_{SDH}','J_{SCAS}');



% % Assume the vectors are already defined:
% % J_SDH: Nx1 vector
% % J_SCAS: Nx1 vector
% % J_Tot: Nx1 vector
% 
% % Step 1: Create grid for surface
% x = SDH(Sm,Mm);   % Ensure column vector
% y = SCAS(Sm,Pm);  % Ensure column vector
% z = J_M(:);   % Ensure column vector
% 
% % Step 2: Define a grid for interpolation
% [xq, yq] = meshgrid(linspace(min(x), max(x), 100), ...
%                     linspace(min(y), max(y), 100));
% 
% % Step 3: Interpolate scattered data onto grid
% zq = griddata(x, y, z, xq, yq, 'cubic');  % or 'linear', 'nearest'
% 
% % Step 4: Plot surface
% figure;
% surf(xq, yq, zq);
% xlabel('J_{SDH} (\mumol/min/mg)');
% ylabel('J_{SCAS} (\mumol/min/mg)');
% zlabel('J_{Tot} (\mumol/min/mg)');
% title('Surface Plot of J_{Tot} as a Function of J_{SDH} and J_{SCAS}');
% shading interp;  % Smooth surface
% colorbar;
% grid on;
% set(gca, 'FontSize', 12);
% 
% [X, Y] = meshgrid(SDH(Sm,Mm),SCAS(Sm,Pm));
% [Z] = meshgrid(J_S);
% figure;
% surf(X, Y, Z);
% xlabel('J_{SDH} (\mumol/min/mg)');
% ylabel('J_{SCAS} (\mumol/min/mg)');
% zlabel('J_{Tot} (\mumol/min/mg)');
% title('Surface Plot of J_{Tot} as a Function of J_{SDH} and J_{SCAS}');
% shading interp;  % Smooth surface
% colorbar;
% grid on;
% set(gca, 'FontSize', 12);
% 
% 
% % Visualize J_S vs J_SDH and J_S vs J_SCAS
% figure('Position', [100, 100, 1200, 400]);
% 
% % Plot J_S vs J_SDH
% subplot(1, 2, 1);
% scatter(SDH(Sm,Mm), J_S, 50, Time, 'filled');
% colormap jet;
% colorbar;
% xlabel('J_{SDH} (µmol/min/g)');
% ylabel('J_S (µmol/min/g)');
% title('Net Succinate Flux vs. SDH Flux');
% grid on;
% 
% % Plot J_S vs J_SCAS
% subplot(1, 2, 2);
% scatter(SCAS(Sm,Pm), J_S, 50, Time, 'filled');
% colormap jet;
% colorbar;
% xlabel('J_{SCAS} (µmol/min/g)');
% ylabel('J_S (µmol/min/g)');
% title('Net Succinate Flux vs. SCAS Flux');
% grid on;
% 
% 
% % 3D Scatter Plot: J_S vs J_SDH and J_SCAS
% figure('Position', [100, 100, 800, 600]);
% scatter3(SDH(Sm,Mm), SCAS(Sm,Pm), J_S, 50, Time, 'filled');
% colormap jet;
% colorbar;
% xlabel('J_{SDH} (µmol/min/g)', 'FontSize', 12);
% ylabel('J_{SCAS} (µmol/min/g)', 'FontSize', 12);
% zlabel('J_S (µmol/min/g)', 'FontSize', 12);
% title('3D Plot of Net Succinate Flux vs. SDH and SCAS Fluxes', 'FontSize', 14);
% grid on;
% view(-45, 30); % Adjust view angle for better visualization
% 
% J_SDH=SDH(Sm,Mm);
% J_SCAS=SCAS(Sm,Pm);
% % Create bins for J_SDH and J_SCAS
% n_bins = 50;
% J_SDH_edges = linspace(min(J_SDH), max(J_SDH), n_bins + 1);
% J_SCAS_edges = linspace(min(J_SCAS), max(J_SCAS), n_bins + 1);
% J_S_grid = zeros(n_bins, n_bins);
% 
% % Bin data and compute average J_S
% for i = 1:n_bins
%     for j = 1:n_bins
%         % Find indices where J_SDH and J_SCAS fall into the current bin
%         idx = J_SDH >= J_SDH_edges(i) & J_SDH < J_SDH_edges(i+1) & ...
%               J_SCAS >= J_SCAS_edges(j) & J_SCAS < J_SCAS_edges(j+1);
%         if any(idx)
%             J_S_grid(i, j) = mean(J_S(idx), 'omitnan');
%         else
%             J_S_grid(i, j) = NaN; % No data in this bin
%         end
%     end
% end
% 
% % Create heatmap
% figure('Position', [100, 100, 800, 600]);
% h = heatmap(J_SDH_edges(1:end-1), J_SCAS_edges(1:end-1), J_S_grid', ...
%     'Colormap', jet, ...
%     'ColorbarVisible', 'on', ...
%     'Title', 'Heatmap of Net Succinate Flux (J_S) vs. SDH and SCAS Fluxes', ...
%     'XLabel', 'J_{SDH} (µmol/min/g)', ...
%     'YLabel', 'J_{SCAS} (µmol/min/g)', ...
%     'GridVisible', 'off');
% % h.Colorbar.Label.String = 'J_S (µmol/min/g)';
% h.Colorbar.Label.FontSize = 12;




% X1=Time;
% Y1=Mm;
% Y2=Mc;
% Y3=Sm;
% Y4=Sc;
% Y5=Pm;
% Y6=Pc;
% createfigure(X1, Y1, Y2, Y3, Y4, Y5, Y6)
% function createfigure(X1, Y1, Y2, Y3, Y4, Y5, Y6)
% %CREATEFIGURE(X1, Y1, Y2, Y3, Y4, Y5, Y6)
% %  X1:  vector of plot x data
% %  Y1:  vector of plot y data
% %  Y2:  vector of plot y data
% %  Y3:  vector of plot y data
% %  Y4:  vector of plot y data
% %  Y5:  vector of plot y data
% %  Y6:  vector of plot y data
% 
% %  Auto-generated by MATLAB on 09-Jun-2025 17:15:50
% 
% % Create figure
% figure('InvertHardcopy','off','Color',[1 1 1]);
% 
% % Create subplot
% subplot1 = subplot(2,3,1);
% hold(subplot1,'on');
% 
% % Create plot
% plot(X1,Y1,'DisplayName','M_m (mitochondria)','LineWidth',3,'Color',[0 0 0]);
% 
% % Create title
% title('[M_m]','FontSize',33);
% 
% % Uncomment the following line to preserve the X-limits of the axes
% % xlim(subplot1,[0 1]);
% % Uncomment the following line to preserve the Y-limits of the axes
% % ylim(subplot1,[10 12.4331321136795]);
% % Uncomment the following line to preserve the Z-limits of the axes
% % zlim(subplot1,[-1 1]);
% % box(subplot1,'on');
% grid(subplot1,'on');
% hold(subplot1,'off');
% % Set the remaining axes properties
% set(subplot1,'FontSize',30);
% % Create subplot
% subplot2 = subplot(2,3,4);
% hold(subplot2,'on');
% 
% % Create plot
% plot(X1,Y2,'DisplayName','M_m (mitochondria)','LineWidth',3,'Color',[0 0 0]);
% 
% % Create title
% title('[M_c]','FontSize',33);
% 
% % Uncomment the following line to preserve the X-limits of the axes
% % xlim(subplot2,[0 1]);
% % Uncomment the following line to preserve the Y-limits of the axes
% % ylim(subplot2,[2 4.8688970607931]);
% % Uncomment the following line to preserve the Z-limits of the axes
% % zlim(subplot2,[-1 1]);
% % box(subplot2,'on');
% grid(subplot2,'on');
% hold(subplot2,'off');
% % Set the remaining axes properties
% set(subplot2,'FontSize',30);
% % Create subplot
% subplot3 = subplot(2,3,2);
% hold(subplot3,'on');
% 
% % Create plot
% plot(X1,Y3,'DisplayName','S_m (mitochondria)','LineWidth',3,'Color',[0 0 0]);
% 
% % Create title
% title('[S_m]','FontSize',33);
% 
% % Uncomment the following line to preserve the X-limits of the axes
% % xlim(subplot3,[0 1]);
% % Uncomment the following line to preserve the Y-limits of the axes
% % ylim(subplot3,[0 10]);
% % Uncomment the following line to preserve the Z-limits of the axes
% % zlim(subplot3,[-1 1]);
% % box(subplot3,'on');
% grid(subplot3,'on');
% hold(subplot3,'off');
% % Set the remaining axes properties
% set(subplot3,'FontSize',30);
% % Create subplot
% subplot4 = subplot(2,3,5);
% hold(subplot4,'on');
% 
% % Create plot
% plot(X1,Y4,'DisplayName','S_m (mitochondria)','LineWidth',3,'Color',[0 0 0]);
% 
% % Create title
% title('[S_c]','FontSize',33);
% 
% % Uncomment the following line to preserve the X-limits of the axes
% % xlim(subplot4,[0 1]);
% % Uncomment the following line to preserve the Y-limits of the axes
% % ylim(subplot4,[0.323869964341073 1.0834655421181]);
% % Uncomment the following line to preserve the Z-limits of the axes
% % zlim(subplot4,[-1 1]);
% % box(subplot4,'on');
% grid(subplot4,'on');
% hold(subplot4,'off');
% % Set the remaining axes properties
% set(subplot4,'FontSize',30);
% % Create subplot
% subplot5 = subplot(2,3,3);
% hold(subplot5,'on');
% 
% % Create plot
% plot(X1,Y5,'DisplayName','P_m (mitochondria)','LineWidth',3,'Color',[0 0 0]);
% 
% % Create title
% title('[P_m]','FontSize',33);
% 
% % Uncomment the following line to preserve the X-limits of the axes
% % xlim(subplot5,[0 1]);
% % Uncomment the following line to preserve the Y-limits of the axes
% % ylim(subplot5,[0 10]);
% % Uncomment the following line to preserve the Z-limits of the axes
% % zlim(subplot5,[-1 1]);
% % box(subplot5,'on');
% grid(subplot5,'on');
% hold(subplot5,'off');
% % Set the remaining axes properties
% set(subplot5,'FontSize',30);
% % Create subplot
% subplot6 = subplot(2,3,6);
% hold(subplot6,'on');
% 
% % Create plot
% plot(X1,Y6,'DisplayName','P_m (mitochondria)','LineWidth',3,'Color',[0 0 0]);
% 
% % Create title
% title('[P_c]','FontSize',33);
% 
% % Uncomment the following line to preserve the X-limits of the axes
% % xlim(subplot6,[0 1]);
% % Uncomment the following line to preserve the Y-limits of the axes
% % ylim(subplot6,[0.5 2.36134511569117]);
% % Uncomment the following line to preserve the Z-limits of the axes
% % zlim(subplot6,[-1 1]);
% % box(subplot6,'on');
% grid(subplot6,'on');
% hold(subplot6,'off');
% % Set the remaining axes properties
% set(subplot6,'FontSize',30);
% end