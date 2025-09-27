%% =========================
% Prior vs Posterior Plots
% =========================

load('MCMC_Result.mat')

% Number of prior samples for plotting
n_prior = 100000;

% --- Generate prior samples ---
% First 5 parameters: truncated Gaussian >0
prior_Ts = draw_pos_norm(mu_vec(1), sigma_vec(1), n_prior);
prior_Tm = draw_pos_norm(mu_vec(2), sigma_vec(2), n_prior);
prior_Ks = draw_pos_norm(mu_vec(3), sigma_vec(3), n_prior);
prior_Km = draw_pos_norm(mu_vec(4), sigma_vec(4), n_prior);
prior_Kp = draw_pos_norm(mu_vec(5), sigma_vec(5), n_prior);

% Lambdas: uniform(a_c,b_c)
prior_lambda21 = unifrnd(a_c, b_c, n_prior, 1);
prior_lambda31 = unifrnd(a_c, b_c, n_prior, 1);

% Sigma^2: Inverse-Gamma(alpha0,beta0)
prior_sigma2 = 1 ./ gamrnd(alpha_0, 1/beta_0, n_prior, 1);

% --- Posterior draws ---
post_Ts       = theta_samples(:,1);
post_Tm       = theta_samples(:,2);
post_Ks       = theta_samples(:,3);
post_Km       = theta_samples(:,4);
post_Kp       = theta_samples(:,5);
post_lambda21 = theta_samples(:,6);
post_lambda31 = theta_samples(:,7);
post_sigma2   = theta_samples(:,8);

% --- Plotting ---
% param_names = { 'T^{s}_{\\max}','T^{m}_{\\max}', ...
%     'K^{s}_{m}','K^{m}_{m}','K^{p}_{m}', ...
%     '\\lambda_{21}','\\lambda_{31}','\\sigma^2'};
param_names = { ...
    'T^s_{max}', ...
    'T^m_{max}', ...
    'K^s_m', ...
    'K^m_m', ...
    'K^p_m', ...
    '\lambda_{21}', ...
    '\lambda_{31}','\sigma^2'};

prior_samples = {prior_Ts, prior_Tm, prior_Ks, prior_Km, prior_Kp,...
                 prior_lambda21, prior_lambda31, prior_sigma2};
post_samples  = {post_Ts, post_Tm, post_Ks, post_Km, post_Kp,...
                 post_lambda21, post_lambda31, post_sigma2};

figure;
for k = 1:8
    subplot(2,4,k); hold on;
    histogram(prior_samples{k},50,'Normalization','pdf','FaceColor',[0.3 0.3 1],'EdgeColor','none','FaceAlpha',0.5,'DisplayName','Prior');
    histogram(post_samples{k},50,'Normalization','pdf','FaceColor',[1 0.3 0.3],'EdgeColor','none','FaceAlpha',0.5,'DisplayName','Posterior');
    xlabel(param_names{k},'Interpreter','tex','FontSize',11);
    ylabel('Density');
    if k == 1
        legend('show');
    end
    % title(['Prior vs Posterior: ',param_names{k}],'Interpreter','latex');
    hold off;
end

%% --- Helper function for truncated Gaussian draws (>0) ---
function X = draw_pos_norm(mu, sd, n)
    X = mu + sd.*randn(n,1);
    while any(X<=0)
        idx = (X<=0);
        X(idx) = mu + sd.*randn(sum(idx),1);
    end
end
