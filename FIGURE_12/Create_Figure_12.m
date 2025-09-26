%% 3x3 figure: concentrations (m & ims) + fluxes (no J_tot)
% Expects both MAT files to contain:
% Time, Mm, Mims, Sm, Sims, Pm, Pims, J_M, J_S, J_P

clear; clc;

% -------- Load data --------
SCR = load('SDH_SCR_Simulated_Data.mat');  % SDH Active
KO  = load('SDH_KO_simulated_Data.mat');   % SDH Deficient

t1 = SCR.Time;  t2 = KO.Time;

% -------- Figure & layout --------
fig = figure('Color',[1 1 1], 'InvertHardcopy','off', 'Name','3x3: States + Fluxes (no J_{tot})');
tl  = tiledlayout(fig, 3, 3, 'TileSpacing','compact', 'Padding','compact');

% Helpers
lw = 2.3; fs = 13;
add_label = @(ax,lab) text(ax, 0.02, 0.93, lab, 'Units','normalized', ...
    'FontWeight','bold','FontSize',12);

% ===== Row 1: matrix concentrations =====
% (a) [M_m]
ax1 = nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on');
plot(ax1, t1, SCR.Mm, 'LineWidth', lw, 'DisplayName','SDH Active');
plot(ax1, t2, KO.Mm,  'LineWidth', lw, 'LineStyle','--', 'DisplayName','SDH Deficient');
xlabel(ax1,'Time (min)','FontWeight','bold'); ylabel(ax1,'[M_m] (mM)','FontWeight','bold');
set(ax1,'FontSize',fs,'FontWeight','bold'); add_label(ax1,'a)');

% (b) [S_m]
ax2 = nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on');
plot(ax2, t1, SCR.Sm, 'LineWidth', lw);
plot(ax2, t2, KO.Sm,  'LineWidth', lw, 'LineStyle','--');
xlabel(ax2,'Time (min)','FontWeight','bold'); ylabel(ax2,'[S_m] (mM)','FontWeight','bold');
set(ax2,'FontSize',fs,'FontWeight','bold'); add_label(ax2,'b)');

% (c) [P_m]
ax3 = nexttile(tl,3); hold(ax3,'on'); grid(ax3,'on');
plot(ax3, t1, SCR.Pm, 'LineWidth', lw);
plot(ax3, t2, KO.Pm,  'LineWidth', lw, 'LineStyle','--');
xlabel(ax3,'Time (min)','FontWeight','bold'); ylabel(ax3,'[P_m] (mM)','FontWeight','bold');
set(ax3,'FontSize',fs,'FontWeight','bold'); add_label(ax3,'c)');

% ===== Row 2: IMS concentrations =====
% (d) [M_{ims}]
ax4 = nexttile(tl,4); hold(ax4,'on'); grid(ax4,'on');
plot(ax4, t1, SCR.Mims, 'LineWidth', lw);
plot(ax4, t2, KO.Mims,  'LineWidth', lw, 'LineStyle','--');
xlabel(ax4,'Time (min)','FontWeight','bold'); ylabel(ax4,'[M_{ims}] (mM)','FontWeight','bold');
set(ax4,'FontSize',fs,'FontWeight','bold'); add_label(ax4,'d)');

% (e) [S_{ims}]
ax5 = nexttile(tl,5); hold(ax5,'on'); grid(ax5,'on');
plot(ax5, t1, SCR.Sims, 'LineWidth', lw);
plot(ax5, t2, KO.Sims,  'LineWidth', lw, 'LineStyle','--');
xlabel(ax5,'Time (min)','FontWeight','bold'); ylabel(ax5,'[S_{ims}] (mM)','FontWeight','bold');
set(ax5,'FontSize',fs,'FontWeight','bold'); add_label(ax5,'e)');

% (f) [P_{ims}]
ax6 = nexttile(tl,6); hold(ax6,'on'); grid(ax6,'on');
plot(ax6, t1, SCR.Pims, 'LineWidth', lw);
plot(ax6, t2, KO.Pims,  'LineWidth', lw, 'LineStyle','--');
xlabel(ax6,'Time (min)','FontWeight','bold'); ylabel(ax6,'[P_{ims}] (mM)','FontWeight','bold');
set(ax6,'FontSize',fs,'FontWeight','bold'); add_label(ax6,'f)');

% ===== Row 3: Fluxes (no total) =====
% (g) J_mal
ax7 = nexttile(tl,7); hold(ax7,'on'); grid(ax7,'on');
plot(ax7, t1, SCR.J_M, 'LineWidth', lw);
plot(ax7, t2, KO.J_M,  'LineWidth', lw, 'LineStyle','--');
xlabel(ax7,'Time (min)','FontWeight','bold'); ylabel(ax7,{'J_{mal}','(\mumol/min/g)'},'FontWeight','bold');
set(ax7,'FontSize',fs,'FontWeight','bold'); add_label(ax7,'g)');

% (h) J_succ
ax8 = nexttile(tl,8); hold(ax8,'on'); grid(ax8,'on');
plot(ax8, t1, SCR.J_S, 'LineWidth', lw);
plot(ax8, t2, KO.J_S,  'LineWidth', lw, 'LineStyle','--');
xlabel(ax8,'Time (min)','FontWeight','bold'); ylabel(ax8,{'J_{succ}','(\mumol/min/g)'},'FontWeight','bold');
set(ax8,'FontSize',fs,'FontWeight','bold'); add_label(ax8,'h)');

% (i) J_pho
ax9 = nexttile(tl,9); hold(ax9,'on'); grid(ax9,'on');
plot(ax9, t1, SCR.J_P, 'LineWidth', lw);
plot(ax9, t2, KO.J_P,  'LineWidth', lw, 'LineStyle','--');
xlabel(ax9,'Time (min)','FontWeight','bold'); ylabel(ax9,{'J_{pho}','(\mumol/min/g)'},'FontWeight','bold');
set(ax9,'FontSize',fs,'FontWeight','bold'); add_label(ax9,'i)');

% % -------- Shared legend (top, horizontal) --------
% lg = legend(tl, {'SDH Active','SDH Deficient'});
% lg.Layout.Tile = 'north';
% lg.Orientation  = 'horizontal';
% lg.FontSize = 11;
% lg.EdgeColor = [1 1 1];

% -------- (Optional) axis limits per panel --------
% xlim(ax1,[0 2.5]); xlim(ax2,[0 2.5]); xlim(ax3,[0 2.5]);
% xlim(ax4,[0 2.5]); xlim(ax5,[0 2.5]); xlim(ax6,[0 2.5]);
% xlim(ax7,[0 2.5]); xlim(ax8,[0 2.5]); xlim(ax9,[0 2.5]);
% ylim(ax7,[-2 100]);   % J_mal
% ylim(ax8,[-10 2]);    % J_succ
% ylim(ax9,[-75 5]);    % J_pho

% -------- (Optional) save --------
% exportgraphics(fig, 'figure_3x3_states_fluxes_no_Jtot.png', 'Resolution', 300);
