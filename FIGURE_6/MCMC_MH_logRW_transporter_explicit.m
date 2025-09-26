% function MCMC_MH_logRW_transporter_explicit()
% ============================================================
% Metropolis–Hastings MCMC for the nonlinear transporter model
%  - Log-normal random-walk proposals (positivity preserved)
%  - Explicit log-likelihood + log-prior in log_alpha
%  - Correct Hastings correction for asymmetry
%  - Mixed priors:
%       * Trunc-Gaussian (>0) for first 5 params
%       * Uniform for lambda21, lambda31
%       * Inverse-Gamma for sigma^2 (obs variance)
%  - Uses your obj_func_combined() flux equations (SSE likelihood)
%  - Posterior flux bands + prior-vs-posterior + traces
% ============================================================
clear; clc; close all;
rng(42);

%% ----------------- Step 1: Load Experimental Data -----------------
% Succinate data
load('Palmier_DIC_Exp_Data.mat')
Pm_v_succ = 10 * ones(1, 6);
Sm_v_succ = 10 * ones(1, 6);
Mm_v_succ = 10 * ones(1, 6);
Pc_v_succ = 0 * ones(1, 6);
Mc_v_succ = 0 * ones(1, 6);
Sc_v_succ = S;
J_succ_exp = T_S; % Experimental data for J_succ

% Malate data
load('Palmier_DIC_Competition_Exp_Data.mat')
Pm_v_mal = 10 * ones(1, 6);
Sm_v_mal = 10 * ones(1, 6);
Mm_v_mal = 10 * ones(1, 6);
Pc_v_mal = 0 * ones(1, 6);
Mc_v_mal = 1 ./ M_with_S;
Sc_v_mal = 0.5 * ones(1, 6);
J_mal_exp = 1 ./ T_M_with_S; % Experimental data for J_mal

% Phosphate data
load('Palmier_DIC_Pho_Exp_Data.mat')
Pm_v_pho = 10 * ones(1, 6);
Sm_v_pho = 10 * ones(1, 6);
Mm_v_pho = 10 * ones(1, 6);
Pc_v_pho = 0 * ones(1, 6);
Mc_v_pho =  1 ./ M_control;
Sc_v_pho = 0 * ones(1, 6);
J_pho_exp = 1 ./ T_M_control; % Experimental data for J_pho

% % Phosphate data
% load('Palmier_DIC_Pho_Exp_Data.mat')
% Pm_v_pho = 10 * ones(1, 6);
% Sm_v_pho = 10 * ones(1, 6);
% Mm_v_pho = 10 * ones(1, 6);
% Pc_v_pho = 0 * ones(1, 6);
% Mc_v_pho =  P_with_M(1,1:6);
% Sc_v_pho = 0 * ones(1, 6);
% J_pho_exp = T_P_with_M(1,1:6); % Experimental data for J_pho

% Volumes (assumed)
Vm   = 1;
Vims = 1/10;

% Pack data into a struct for easy passing
data = struct( ...
   'Vims',Vims,'Vm',Vm, ...
   'Mc_v_succ',Mc_v_succ,'Sm_v_succ',Sm_v_succ,'Sc_v_succ',Sc_v_succ,'Mm_v_succ',Mm_v_succ,'Pm_v_succ',Pm_v_succ,'Pc_v_succ',Pc_v_succ,'J_succ_exp',J_succ_exp, ...
   'Mc_v_mal',Mc_v_mal,'Sm_v_mal',Sm_v_mal,'Sc_v_mal',Sc_v_mal,'Mm_v_mal',Mm_v_mal,'Pm_v_mal',Pm_v_mal,'Pc_v_mal',Pc_v_mal,'J_mal_exp',J_mal_exp, ...
   'Mc_v_pho',Mc_v_pho,'Sm_v_pho',Sm_v_pho,'Sc_v_pho',Sc_v_pho,'Mm_v_pho',Mm_v_pho,'Pm_v_pho',Pm_v_pho,'Pc_v_pho',Pc_v_pho,'J_pho_exp',J_pho_exp);

n_total = 18; % 6 succinate + 6 malate + 6 phosphate

%% ----------------- Step 2: Priors (your settings) -----------------
% First five parameters (positive): [Ts_max, Tm_max, Ks_m, Km_m, Kp_m]
mu_vec    = [64.2;  69.04;  1.17; 0.23; 0.93];   % means
sigma_vec = [5.1;    5.6;  0.10; 0.02; 0.10];   % SDs

% Uniform priors for last two parameters (lambda21, lambda31)
a_c = 0.01; b_c = 50;

% Inverse-Gamma prior for sigma^2
alpha_0 = 1; beta_0 = 1;

% Prebuild prior/likelihood handles
log_prior_all = @(params, sigma2) log_prior_mix(params, sigma2, mu_vec, sigma_vec, a_c, b_c, alpha_0, beta_0);
log_likelihood = @(params, sigma2) -n_total/2*log(2*pi*sigma2) - 0.5*(1/sigma2)* obj_sse(params, data);

%% ----------------- Step 3: Metropolis-Hastings --------------------
% Parameter vector:
%   [Ts_max, Tm_max, Ks_m, Km_m, Kp_m, lambda21, lambda31, sigma2]
theta = [59; 61; 1.0; 0.3; 1.41; 1; 0.51; 0.23];  % initial values (>0)

% MCMC controls
n_samples = 100000;      % total iterations
burn_in   = 2000;       % discard first
thin      = 1;          % keep every thin-th sample
z_step    = 0.01;        % log-RW step size (std of log increment)
report_every = 2000;

p_all = numel(theta);
theta_samples_raw = zeros(n_samples, p_all);
acc = false(n_samples,1);

% ---- Cache current log-likelihood and log-prior explicitly ----
loglik_curr   = log_likelihood(theta(1:7), theta(8));
logprior_curr = log_prior_all(theta(1:7), theta(8));

for i = 1:n_samples
    % Log-normal random walk proposal for all 8 positive entries
    log_theta_prop = log(theta) + z_step * randn(p_all,1);
    theta_prop = exp(log_theta_prop);

    % Reject if uniforms out of range or sigma2 <= 0 (exp keeps it >0)
    if any(theta_prop(6:7) < a_c | theta_prop(6:7) > b_c) || theta_prop(8) <= 0
        theta_samples_raw(i,:) = theta.'; 
        continue;
    end

    % Proposed log-likelihood and log-prior
    loglik_prop    = log_likelihood(theta_prop(1:7), theta_prop(8));
    logprior_prop  = log_prior_all(theta_prop(1:7), theta_prop(8));

    % Hastings correction for asymmetric log-normal proposal:
    % log q(curr|prop) - log q(prop|curr) = sum(log(theta) - log(theta_prop))
    log_hastings = sum(log(theta) - log(theta_prop));

    % MH acceptance on log-scale (explicit terms)
    log_alpha = (loglik_prop - loglik_curr) + (logprior_prop - logprior_curr) + log_hastings;

    if log(rand) < log_alpha
        theta         = theta_prop;      % accept
        loglik_curr   = loglik_prop;     % cache for next iteration
        logprior_curr = logprior_prop;
        acc(i)        = true;
    end

    theta_samples_raw(i,:) = theta.';
    if mod(i, report_every)==0
        fprintf('Iter %7d | acc(last %d)=%.3f | acc(all)=%.3f | sigma^2=%.4g\n', ...
            i, report_every, mean(acc(i-report_every+1:i)), mean(acc(1:i)), theta(8));
    end
end
fprintf('\nOverall acceptance rate: %.3f\n', mean(acc));

% Burn-in + thinning
keep_idx = (burn_in+1):thin:n_samples;
theta_samples = theta_samples_raw(keep_idx,:);

%% ----------------- Step 4: Posterior Flux Predictions --------------
n_post = size(theta_samples, 1);
J_S_pred = zeros(n_post, 6);
J_M_pred = zeros(n_post, 6);
J_P_pred = zeros(n_post, 6);

for i = 1:n_post
    params = theta_samples(i, 1:7);
    fluxes = forward_fluxes(params, data); % returns {J_S, J_M, J_P}
    J_S_pred(i, :) = fluxes{1};
    J_M_pred(i, :) = fluxes{2};
    J_P_pred(i, :) = fluxes{3};
end

% Posterior means & 95% CIs
mean_J_S = mean(J_S_pred, 1);
mean_J_M = mean(J_M_pred, 1);
mean_J_P = mean(J_P_pred, 1);
CI_J_S = prctile(J_S_pred, [2.5, 97.5], 1);
CI_J_M = prctile(J_M_pred, [2.5, 97.5], 1);
CI_J_P = prctile(J_P_pred, [2.5, 97.5], 1);

%% ----------------- Step 5: Plots - Flux Bands ---------------------
figure('Color','w','Name','Flux credible bands');
% Succinate
subplot(1, 3, 1); hold on;
fill([Sc_v_succ, fliplr(Sc_v_succ)], [CI_J_S(1,:), fliplr(CI_J_S(2,:))], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% CI');
plot(Sc_v_succ, mean_J_S, 'r--', 'LineWidth', 2, 'DisplayName', 'Posterior mean');
plot(Sc_v_succ, J_succ_exp, 'b-', 'LineWidth', 2, 'DisplayName', 'Experimental');
xlabel('Succinate [S_c]'); ylabel('J_S'); title('Succinate Flux'); legend; grid on; hold off;

% Malate
subplot(1, 3, 2); hold on;
fill([Mc_v_mal, fliplr(Mc_v_mal)], [CI_J_M(1,:), fliplr(CI_J_M(2,:))], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% CI');
plot(Mc_v_mal, mean_J_M, 'r--', 'LineWidth', 2, 'DisplayName', 'Posterior mean');
plot(Mc_v_mal, J_mal_exp, 'b-', 'LineWidth', 2, 'DisplayName', 'Experimental');
xlabel('Malate [M_c]'); ylabel('J_M'); title('Malate Flux'); legend; grid on; hold off;

% Phosphate
subplot(1, 3, 3); hold on;
fill([Mc_v_pho, fliplr(Mc_v_pho)], [CI_J_P(1,:), fliplr(CI_J_P(2,:))], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% CI');
plot(Mc_v_pho, mean_J_P, 'r--', 'LineWidth', 2, 'DisplayName', 'Posterior mean');
plot(Mc_v_pho, J_pho_exp, 'b-', 'LineWidth', 2, 'DisplayName', 'Experimental');
xlabel('Phosphate [P_c]'); ylabel('J_P'); title('Phosphate Flux'); legend; grid on; hold off;

% Print CIs
fprintf('95%% Credible Intervals for J_S:\n');
for i = 1:length(Sc_v_succ)
    fprintf('S_c = %.3f: [%.4f, %.4f]\n', Sc_v_succ(i), CI_J_S(1,i), CI_J_S(2,i));
end
fprintf('95%% Credible Intervals for J_M:\n');
for i = 1:length(Mc_v_mal)
    fprintf('M_c = %.3f: [%.4f, %.4f]\n', Mc_v_mal(i), CI_J_M(1,i), CI_J_M(2,i));
end
fprintf('95%% Credible Intervals for J_P:\n');
for i = 1:length(Mc_v_pho)
    fprintf('P_c = %.3f: [%.4f, %.4f]\n', Mc_v_pho(i), CI_J_P(1,i), CI_J_P(2,i));
end

%% ----------------- Step 6: Prior draws (for overlays) --------------
n_prior = 10000;
Tmax_sf_prior  = draw_pos_norm(mu_vec(1), sigma_vec(1), n_prior);
Tmax_mf_prior  = draw_pos_norm(mu_vec(2), sigma_vec(2), n_prior);
km_sf_prior    = draw_pos_norm(mu_vec(3), sigma_vec(3), n_prior);
km_mf_prior    = draw_pos_norm(mu_vec(4), sigma_vec(4), n_prior);
km_pf_prior    = draw_pos_norm(mu_vec(5), sigma_vec(5), n_prior);
lambda21_prior = unifrnd(a_c, b_c, n_prior, 1);
lambda31_prior = unifrnd(a_c, b_c, n_prior, 1);
sigma2_prior   = 1 ./ gamrnd(alpha_0, 1/beta_0, n_prior, 1);

%% ----------------- Step 7: Prior vs Posterior ----------------------
post = theta_samples; % shorthand
labels = {'T^{s}_{max}','T^{m}_{max}','k^{s}_m','k^{m}_m','k^{p}_m','\lambda_{21}','\lambda_{31}','\sigma^2'};
priors = {Tmax_sf_prior, Tmax_mf_prior, km_sf_prior, km_mf_prior, km_pf_prior, lambda21_prior, lambda31_prior, sigma2_prior};

figure('Color','w','Name','Prior vs Posterior');
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');
for k=1:8
    nexttile; hold on;
    histogram(priors{k}, 60, 'Normalization','pdf','FaceColor',[0.2 0.5 1], 'EdgeColor','none','DisplayName','Prior');
    histogram(post(:,k),  60, 'Normalization','pdf','FaceColor',[1 0.3 0.3], 'EdgeColor','none','FaceAlpha',0.7,'DisplayName','Posterior');
    title(labels{k},'Interpreter','tex','FontWeight','bold');
    legend('Location','best'); grid on; hold off;
end

%% ----------------- Step 8: Traces & Posterior summary --------------
figure('Color','w','Name','MCMC traces');
tiledlayout(4,2,'Padding','compact','TileSpacing','compact');
for k=1:8
    nexttile; plot(theta_samples_raw(:,k), 'LineWidth', 0.8);
    xlabel('iteration'); ylabel(labels{k},'Interpreter','tex'); grid on;
end

posterior_mean = mean(theta_samples, 1);
posterior_std  = std(theta_samples, 0, 1);
fprintf(['\nPosterior means:\nTs_max=%.3f, Tm_max=%.3f, Ks_m=%.3f, Km_m=%.3f, Kp_m=%.3f, ' ...
         'lambda21=%.3f, lambda31=%.3f, sigma2=%.4g\n'], ...
         posterior_mean(1), posterior_mean(2), posterior_mean(3), posterior_mean(4), ...
         posterior_mean(5), posterior_mean(6), posterior_mean(7), posterior_mean(8));
fprintf(['Posterior stds :\nTs_max=%.3f, Tm_max=%.3f, Ks_m=%.3f, Km_m=%.3f, Kp_m=%.3f, ' ...
         'lambda21=%.3f, lambda31=%.3f, sigma2=%.4g\n'], ...
         posterior_std(1), posterior_std(2), posterior_std(3), posterior_std(4), ...
         posterior_std(5), posterior_std(6), posterior_std(7), posterior_std(8));
% end
% ========================= END MAIN =========================


% ===================== Helper Functions =====================
function sse = obj_sse(params, data)
% Wrapper: calls your flux model and returns SSE across all 3 datasets
flux_or_sse = obj_func_combined(params, data.Vims, data.Vm, ...
    data.Mc_v_succ, data.Sm_v_succ, data.Sc_v_succ, data.Mm_v_succ, data.Pm_v_succ, data.Pc_v_succ, data.J_succ_exp, ...
    data.Mc_v_mal,  data.Sm_v_mal,  data.Sc_v_mal,  data.Mm_v_mal,  data.Pm_v_mal,  data.Pc_v_mal,  data.J_mal_exp, ...
    data.Mc_v_pho,  data.Sm_v_pho,  data.Sc_v_pho,  data.Mm_v_pho,  data.Pm_v_pho,  data.Pc_v_pho,  data.J_pho_exp);
sse = flux_or_sse; % obj_func_combined returns SSE when *_exp are provided
end

function fluxes = forward_fluxes(params, data)
% Calls your model to get {J_S, J_M, J_P} when *_exp are empty
fluxes = obj_func_combined(params, data.Vims, data.Vm, ...
    data.Mc_v_succ, data.Sm_v_succ, data.Sc_v_succ, data.Mm_v_succ, data.Pm_v_succ, data.Pc_v_succ, [], ...
    data.Mc_v_mal,  data.Sm_v_mal,  data.Sc_v_mal,  data.Mm_v_mal,  data.Pm_v_mal,  data.Pc_v_mal,  [], ...
    data.Mc_v_pho,  data.Sm_v_pho,  data.Sc_v_pho,  data.Mm_v_pho,  data.Pm_v_pho,  data.Pc_v_pho,  []);
end

function lp = log_prior_mix(params, sigma2, mu_vec, sigma_vec, a_c, b_c, alpha_0, beta_0)
% Mixed priors:
%  - Truncated-at-0 Gaussians for first 5 params
%  - Uniform(a_c,b_c) for lambda21 (param 6) and lambda31 (param 7)
%  - Inverse-Gamma(alpha_0, beta_0) for sigma^2
if any(params(1:5) <= 0) || sigma2 <= 0
    lp = -inf; return;
end

% Trunc-N(>0) for first five
lp_truncs = 0;
for i=1:5
    lp_truncs = lp_truncs + log_truncnorm0(params(i), mu_vec(i), sigma_vec(i));
end

% Uniforms for lambda21, lambda31
if params(6) < a_c || params(6) > b_c || params(7) < a_c || params(7) > b_c
    lp = -inf; return;
end

% Inv-Gamma for sigma^2
lp_sig2 = log_inv_gamma(sigma2, alpha_0, beta_0);

lp = lp_truncs + lp_sig2; % uniform adds 0 inside [a_c,b_c]
end

function v = log_truncnorm0(x, mu, sd)
% log N(x|mu,sd) truncated to (0,∞)
if x <= 0, v = -inf; return; end
v = -0.5*log(2*pi*sd^2) - 0.5*((x-mu)/sd)^2 - log(1 - normcdf(0, mu, sd));
end

function v = log_inv_gamma(s2, alpha, beta)
% Inv-Gamma(alpha,beta) density: p(s2) ∝ s2^{-alpha-1} exp(-beta/s2)
if s2 <= 0, v = -inf; return; end
v = (-alpha-1)*log(s2) - beta/s2 - gammaln(alpha) + alpha*log(beta);
end

function X = draw_pos_norm(mu, sd, n)
% Draw from N(mu,sd) truncated at (0,∞) by rejection
X = mu + sd.*randn(n,1);
while any(X<=0)
    idx = (X<=0);
    X(idx) = mu + sd.*randn(sum(idx),1);
end
end

% ---------------- Your original obj_func_combined ----------------
function out = obj_func_combined(params, Vims, Vm, Mc_v_succ, Sm_v_succ, Sc_v_succ, ...
    Mm_v_succ, Pm_v_succ, Pc_v_succ, J_succ_exp, Mc_v_mal, Sm_v_mal, Sc_v_mal, ...
    Mm_v_mal, Pm_v_mal, Pc_v_mal, J_mal_exp, Mc_v_pho, Sm_v_pho, Sc_v_pho, ...
    Mm_v_pho, Pm_v_pho, Pc_v_pho, J_pho_exp)

    % Unpack parameters
    Ts_max   = params(1);  % Succinate Tmax
    Tm_max   = params(2);  % Malate Tmax
    Ks_m     = params(3);  % Succinate Km
    Km_m     = params(4);  % Malate Km
    Kp_m     = params(5);  % Phosphate Km
    lambda21 = params(6);
    lambda31 = params(7);

    % Succinate flux
    J_S = zeros(size(Sm_v_succ));
    for i = 1:length(Sm_v_succ)
        Sm = Sm_v_succ(i);  Sims = Sc_v_succ(i);
        Mm = Mm_v_succ(i);  Mims = Mc_v_succ(i);
        Pm = Pm_v_succ(i);  Pims = Pc_v_succ(i);

        delta1 = 1 + Sm/Ks_m + Mm/Km_m + Pm/Kp_m;
        delta2 = 1 + Sims/Ks_m + Mims/Km_m + Pims/Kp_m;

        phi2_num = Sm + (lambda21)*Mm + ((Vm/Vims)^2)*lambda31*Pm;
        phi2_den = Sims + (lambda21)*Mims + ((Vm/Vims)^2)*lambda31*Pims;
        phi2 = phi2_num / phi2_den;

        J_S(i) = (1/Vm) * (Ts_max * (phi2*Sims - Sm) / (Ks_m * delta1 + Ks_m * phi2 * delta2));
    end

    % Malate flux
    J_M = zeros(size(Mc_v_mal));
    for i = 1:length(Mc_v_mal)
        Sm = Sm_v_mal(i);   Sims = Sc_v_mal(i);
        Mm = Mm_v_mal(i);   Mims = Mc_v_mal(i);
        Pm = Pm_v_mal(i);   Pims = Pc_v_mal(i);

        delta1 = 1 + Sm/Ks_m + Mm/Km_m + Pm/Kp_m;
        delta2 = 1 + Sims/Ks_m + Mims/Km_m + Pims/Kp_m;

        phi2_num = Sm + (lambda21)*Mm + ((Vm/Vims)^2)*lambda31*Pm;
        phi2_den = Sims + (lambda21)*Mims + ((Vm/Vims)^2)*lambda31*Pims;
        phi2 = phi2_num / phi2_den;

        J_M(i) = (1/Vm) * (Tm_max * (phi2*Mims - Mm) / (Km_m * delta1 + Km_m * phi2 * delta2));
    end

    % Phosphate "flux" reusing malate-form
    J_P = zeros(size(Pc_v_pho));
    for i = 1:length(Pc_v_pho)
        Sm = Sm_v_pho(i);   Sims = Sc_v_pho(i);
        Mm = Mm_v_pho(i);   Mims = Mc_v_pho(i);
        Pm = Pm_v_pho(i);   Pims = Pc_v_pho(i);

        delta1 = 1 + Sm/Ks_m + Mm/Km_m + Pm/Kp_m;
        delta2 = 1 + Sims/Ks_m + Mims/Km_m + Pims/Kp_m;

        phi2_num = Sm + (lambda21)*Mm + ((Vm/Vims)^2)*lambda31*Pm;
        phi2_den = Sims + (lambda21)*Mims + ((Vm/Vims)^2)*lambda31*Pims;
        phi2 = phi2_num / phi2_den;
        % phi1_num = Sm + (lambda21)*Mm + (lambda31)*Pm;
        % phi1_den = Sims + (lambda21)*Mims + (lambda31)*Pims;
        % phi1 = phi1_num / phi1_den;

        J_P(i) = (1/Vm) * (Tm_max * (phi2*Mims - Mm) / (Km_m * delta1 + Km_m * phi2 * delta2));
         % J_P(i) = -(1/Vims) * (Tm_max * (Mm - phi1*Mims) / (Km_m * delta1 + Km_m * phi1 * delta2));
    end

    % SSE for likelihood, or return fluxes if no experimental arrays provided
    if ~isempty(J_succ_exp) || ~isempty(J_mal_exp) || ~isempty(J_pho_exp)
        err_s = J_S - J_succ_exp;
        err_m = J_M - J_mal_exp;
        err_p = J_P - J_pho_exp;
        out = sum(err_s.^2 + err_m.^2 + err_p.^2);
    else
        out = {J_S, J_M, J_P};
    end
end
