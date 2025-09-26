run('LSC_Succ.m')
run('LSC_Mal.m')
run('LSC_Pho.m')
% Example usage:
load('LSCs_succ.mat');   % contains mean_LSCs_succ (1x7)
load('LSCs_mal.mat');    % contains mean_LSCs_mal (1x7)
load('LSCs_pho.mat');    % contains mean_LSCs_pho (1x7)

createfigure(mean_LSCs_pho', mean_LSCs_mal', mean_LSCs_succ');