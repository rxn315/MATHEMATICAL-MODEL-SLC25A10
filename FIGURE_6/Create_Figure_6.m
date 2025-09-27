 
load('MCMC_Result.mat')

%% ----------------- Step 5: Plots - Flux Bands (styled like createfigure) ---------------------
% Build patch polygons (x concatenated with its flipped copy; same for y)
X_s = [Sc_v_succ, fliplr(Sc_v_succ)];
Y_s = [CI_J_S(1,:), fliplr(CI_J_S(2,:))];

X_m = [Mc_v_mal,  fliplr(Mc_v_mal)];
Y_m = [CI_J_M(1,:), fliplr(CI_J_M(2,:))];

X_p = [Mc_v_pho,  fliplr(Mc_v_pho)];
Y_p = [CI_J_P(1,:), fliplr(CI_J_P(2,:))];

% Assemble "matrix input to plot" (first col = mean, second col = experimental)
Ymat_s = [mean_J_S(:), J_succ_exp(:)];
Ymat_m = [mean_J_M(:), J_mal_exp(:)];
Ymat_p = [mean_J_P(:), J_pho_exp(:)];

% Create figure (white background; same style)
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);

%% ---- Right panel (c): Phosphate ----
axes1 = axes('Position',[0.757305665349136 0.15241935483871 0.172233201581028 0.686290322580645]);
hold(axes1,'on');

% 95% CI patch (red, 20% alpha)
patch('DisplayName','95% Credible Interval','YData',Y_p,'XData',X_p, ...
    'FaceAlpha',0.2,'FaceColor',[1 0 0],'EdgeColor','none');

% Mean (dashed red) + Experimental (black dots)
h = plot(Mc_v_pho, Ymat_p, 'LineWidth', 2);
set(h(1),'DisplayName','Mean Trajectory','LineStyle','--','Color',[1 0 0]);
set(h(2),'DisplayName','Experimental Data','MarkerFaceColor',[0 0 0], ...
    'MarkerEdgeColor',[0 0 0],'MarkerSize',10,'Marker','o','LineStyle','none','Color',[0 0 1]);

ylabel('Malate Flux (J_{mal})','FontSize',33);
xlabel('[M_{ims}]','FontSize',33);
grid(axes1,'on');  hold(axes1,'off');
set(axes1,'FontSize',30);

%% ---- Middle panel (b): Malate ----
axes2 = axes('Position',[0.463003952569169 0.150806451612903 0.172233201581028 0.694354838709677]);
hold(axes2,'on');

patch('DisplayName','95% Credible Interval','YData',Y_m,'XData',X_m, ...
    'FaceAlpha',0.2,'FaceColor',[1 0 0],'EdgeColor','none');

h = plot(Mc_v_mal, Ymat_m, 'LineWidth', 2);
set(h(1),'DisplayName','Mean Trajectory','LineStyle','--','Color',[1 0 0]);
set(h(2),'DisplayName','Experimental Data','MarkerFaceColor',[0 0 0], ...
    'MarkerEdgeColor',[0 0 0],'MarkerSize',10,'Marker','o','LineStyle','none','Color',[0 0 1]);

ylabel('Malate Flux (J_{mal})','FontSize',33);
xlabel('[M_{ims}]','FontSize',33);
grid(axes2,'on');  hold(axes2,'off');
set(axes2,'FontSize',30);

%% ---- Left panel (a): Succinate ----
axes3 = subplot(1,3,1); % keep as subplot to use legend position trick from createfigure
set(axes3,'Position',[0.12 0.15 0.25 0.70]); % similar spacing
hold(axes3,'on');

patch('DisplayName','95% Credible Interval','YData',Y_s,'XData',X_s, ...
    'FaceAlpha',0.2,'FaceColor',[1 0 0],'EdgeColor','none');

h = plot(Sc_v_succ, Ymat_s, 'LineWidth', 2);
set(h(1),'DisplayName','Mean Trajectory','LineStyle','--','Color',[1 0 0]);
set(h(2),'DisplayName','Experimental Data','MarkerFaceColor',[0 0 0], ...
    'MarkerEdgeColor',[0 0 0],'MarkerSize',10,'Marker','o','LineStyle','none','Color',[0 0 1]);

ylabel('Succinate Flux (J_{succ})','FontSize',33);
xlabel('[S_{ims}]','FontSize',33);
grid(axes3,'on');  hold(axes3,'off');
set(axes3,'FontSize',30);

%% ---- Top legend (horizontal, large, edge = white) ----
legend1 = legend(axes3,'show');
set(legend1,'Position',[0.0351 0.9482 1.0510 0.0583], ... % spans top width
    'Orientation','horizontal','FontSize',27,'EdgeColor',[1 1 1]);

%% ---- Panel letters a), b), c) ----
annotation(figure1,'textbox',[0.0422 0.8672 0.0841 0.0754], ...
    'String','a)','LineWidth',2,'FontWeight','bold','FontSize',30,'EdgeColor','none');
annotation(figure1,'textbox',[0.3723 0.8672 0.0854 0.0754], ...
    'String','b)','LineWidth',2,'FontWeight','bold','FontSize',30,'EdgeColor','none');
annotation(figure1,'textbox',[0.6707 0.8640 0.0841 0.0754], ...
    'String','c)','LineWidth',2,'FontWeight','bold','FontSize',30,'EdgeColor','none');

%% ---- (Optional) keep limits similar to your screenshot ----
% xlim(axes3,[0 1]);   ylim(axes3,[8 35]);    % succinate
% xlim(axes2,[0 1]);   ylim(axes2,[20 50]);   % malate
% xlim(axes1,[0 1]);   ylim(axes1,[15 40]);   % phosphate
