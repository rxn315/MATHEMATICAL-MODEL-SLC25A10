%% ============================================
%  Posterior summary + LaTeX table + lambda plots
%  Mode computed from fitted PDF (Normal/Lognormal/Gamma)
%  Requires: theta_samples (N x 8) post burn-in
%  Optional: accepted_count, total_steps
% =============================================
load('MCMC_Result.mat')

clearvars -except theta_samples accepted_count total_steps
clc

% ---- Sanity checks
if ~exist('theta_samples','var')
    error('theta_samples (N x 8) not found in workspace.');
end
if size(theta_samples,2) ~= 8
    error('theta_samples must have 8 columns: [Ts_max, Tm_max, Ks_m, Km_m, Kp_m, lambda21, lambda31, sigma2].');
end

% ---- Hyperparameter (your RW step size on log-scale)
a_hyper = 0.01;

% ---- Names
param_names_tex = { ...
    'T^{s}_{\\max}','T^{m}_{\\max}', ...
    'K^{s}_{m}','K^{m}_{m}','K^{p}_{m}', ...
    '\\lambda_{21}','\\lambda_{31}','\\sigma^2'};

param_names_axis = { ...
    'T^{s}_{\max}','T^{m}_{\max}', ...
    'K^{s}_{m}','K^{m}_{m}','K^{p}_{m}', ...
    '\lambda_{21}','\lambda_{31}','\sigma^2'};

% ---- Summary arrays
n_params = 8;
means  = zeros(1,n_params);
medians= zeros(1,n_params);
modes  = zeros(1,n_params);          % <-- PDF-based mode
stds   = zeros(1,n_params);
ci_low = zeros(1,n_params);
ci_hi  = zeros(1,n_params);

% keep the fitted PDFs (for optional plotting)
pdfFits(n_params) = struct('pd',[],'name','','support',[]);

for j = 1:n_params
    v = theta_samples(:,j);
    v = v(~isnan(v) & isfinite(v));         % clean

    means(j)    = mean(v);
    medians(j)  = median(v);
    stds(j)     = std(v,0,1);
    prc         = prctile(v,[2.5 97.5]);
    ci_low(j)   = prc(1);
    ci_hi(j)    = prc(2);

    % === MODE from parametric PDF ===
    [pdBest, distName] = fit_best_pdf(v);
    modes(j) = mode_from_pdf(pdBest, v);     % analytic when available
    pdfFits(j).pd = pdBest;
    pdfFits(j).name = distName;
    pdfFits(j).support = [min(v) max(v)];
end

% ---- Acceptance stats (optional)
if ~exist('accepted_count','var') || isempty(accepted_count)
    acc_str = '--';
else
    acc_str = sprintf('%d', accepted_count);
end
if ~exist('total_steps','var') || isempty(total_steps)
    total_steps = size(theta_samples,1);
end
acc_rate = NaN;
if exist('accepted_count','var') && ~isempty(accepted_count) && total_steps > 0
    acc_rate = accepted_count / total_steps;
end

% ---- LaTeX table
fmt = @(x) sprintf('%.5g', x);

header = sprintf([ ...
    '\\begin{table}[ht]\n\\centering\n' ...
    '\\caption{MCMC performance and posterior summary for all parameters (hyperparameter $a=%.3g$).}\n' ...
    '\\label{tab:post_all}\n' ...
    '\\small\n' ...
    '\\begin{tabular}{l c c c c c c c c c}\n\\hline\n' ...
    '\\textbf{Parameter} & Mean & Median & Mode (PDF) & Std & 2.5\\%% & 97.5\\%% & $a$ & Accepted & Total \\\\\n\\hline\n'], a_hyper);

rows = '';
for j = 1:n_params
    rows = [rows, sprintf('%s & %s & %s & %s & %s & %s & %s & %.3g & %s & %d \\\\\n', ...
        param_names_tex{j}, ...
        fmt(means(j)), fmt(medians(j)), fmt(modes(j)), fmt(stds(j)), ...
        fmt(ci_low(j)), fmt(ci_hi(j)), a_hyper, acc_str, total_steps)]; %#ok<AGROW>
end
footer = '\\hline\n\\end{tabular}\n\\end{table}\n';

latex_table = [header, rows, footer];

% ---- Print and save LaTeX
fprintf('%s\n', latex_table);
fid = fopen('posterior_summary_table.tex','w'); fprintf(fid, '%s', latex_table); fclose(fid);

% ---- Plots for lambda21 and lambda31 only: trace + histogram + fitted PDF
ix_l21 = 6; ix_l31 = 7;

figure('Color','w','Position',[100 100 1100 450]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% Trace: lambda21
nexttile; plot(theta_samples(:,ix_l21),'k-');
xlabel('Accepted MCMC steps'); ylabel('\lambda_{21}');
title('\lambda_{21}: trace'); grid on;

% Trace: lambda31
nexttile; plot(theta_samples(:,ix_l31),'k-');
xlabel('Accepted MCMC steps'); ylabel('\lambda_{31}');
title('\lambda_{31}: trace'); grid on;

% Posterior: lambda21
nexttile; hold on;
v = theta_samples(:,ix_l21);
histogram(v,50,'Normalization','pdf','FaceColor',[0.85 0.85 0.85],...
             'EdgeColor','none','DisplayName','Posterior (hist)');
[pd21, name21] = fit_best_pdf(v);
[xi21, f21] = pdf_curve(pd21, v);
plot(xi21, f21, 'r-', 'LineWidth',1.8, 'DisplayName', ['PDF: ' name21]);
mu = mean(v); xline(mu,'--b','LineWidth',1.5,'DisplayName','Mean');
xline(mode_from_pdf(pd21,v),'-k','LineWidth',1.8,'DisplayName','Mode (PDF)');
xlabel('\lambda_{21}'); ylabel('pdf');
title('\lambda_{21}: posterior'); legend('Location','best'); grid on; hold off;

% Posterior: lambda31
nexttile; hold on;
v = theta_samples(:,ix_l31);
histogram(v,50,'Normalization','pdf','FaceColor',[0.85 0.85 0.85],...
             'EdgeColor','none','DisplayName','Posterior (hist)');
[pd31, name31] = fit_best_pdf(v);
[xi31, f31] = pdf_curve(pd31, v);
plot(xi31, f31, 'r-', 'LineWidth',1.8, 'DisplayName', ['PDF: ' name31]);
xline(mode_from_pdf(pd31,v),'-k','LineWidth',1.8,'DisplayName','Mode (PDF)');
mu = mean(v); xline(mu,':','Color',[0.3 0.3 0.3],'LineWidth',1.5,'DisplayName','Mean');
xlabel('\lambda_{31}'); ylabel('pdf');
title('\lambda_{31}: posterior'); legend('Location','best'); grid on; hold off;

% ---- (Optional) print acceptance rate
if ~isnan(acc_rate)
    fprintf('Acceptance: %d / %d (%.1f%%)\n', accepted_count, total_steps, 100*acc_rate);
else
    fprintf('Acceptance: -- / %d (rate unknown)\n', total_steps);
end

%% ===========================
% Helpers: best-fit PDF + mode + curve
% ===========================
function [pdBest, nameBest] = fit_best_pdf(v)
% Try Normal / Lognormal / Gamma (respecting support), pick max log-likelihood
    cand = {}; 
    if all(v > 0)
        cand = {'lognormal','gamma','normal'};
    else
        cand = {'normal'};
    end
    llBest = -inf; pdBest = []; nameBest = '';
    for k = 1:numel(cand)
        try
            pd = fitdist(v, cand{k});
            % log-likelihood on the fitted PDF
            ll = sum(log(pdf(pd, v) + realmin));
            if ll > llBest
                llBest = ll; pdBest = pd; nameBest = pd.DistributionName;
            end
        catch
            % skip failures silently
        end
    end
    % Fallback to Normal if nothing worked
    if isempty(pdBest)
        pdBest = fitdist(v,'normal'); nameBest = pdBest.DistributionName;
    end
end

function m = mode_from_pdf(pd, v)
% Closed-form mode when available; else grid search on pdf(pd,x)
    switch lower(pd.DistributionName)
        case 'normal'
            m = pd.mu;                                    % mode = mean
        case 'lognormal'
            % MATLAB's Lognormal: parameters are mu (log-mean), sigma
            m = exp(pd.mu - pd.sigma^2);
        case 'gamma'
            % MATLAB's Gamma: a=shape, b=scale; mode = (a-1)*b for a>1 (else 0)
            if pd.a > 1
                m = (pd.a - 1) * pd.b;
            else
                m = 0;
            end
        otherwise
            % numeric maximization on a grid covering data support
            [xg, fg] = pdf_curve(pd, v);
            [~,ix] = max(fg); m = xg(ix);
    end
end

function [xi, f] = pdf_curve(pd, v)
% Build a smooth pdf curve over a padded support of the data
    lo = min(v); hi = max(v);
    span = hi - lo; pad = 0.1*max(span, eps);
    if strcmpi(pd.DistributionName,'lognormal') || strcmpi(pd.DistributionName,'gamma')
        % positive support
        lo = max(0, lo - pad); 
        xi = linspace(lo, hi + pad, 600);
    else
        xi = linspace(lo - pad, hi + pad, 600);
    end
    f = pdf(pd, xi);
end
