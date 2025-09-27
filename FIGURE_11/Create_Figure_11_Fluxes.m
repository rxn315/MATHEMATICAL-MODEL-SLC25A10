% Build inputs
X1       = Time(:);
YMatrix1 = [J_S(:), J_M(:), J_P(:)];               % [J_succ, J_mal, J_pho]

xpatch   = [Time(1) Time(end) Time(end) Time(1)];
ymin     = -210;             % choose to cover your data
ymax     =  210;
yl_eff   = [ymin ymin 0 0];  % Efflux region (blue)
yl_inf   = [0 0 ymax ymax];  % Influx region (red)

% Draw the figure
createfigure_Flux(yl_eff, xpatch, yl_inf, X1, YMatrix1);

% (Optional) match the zooms shown in your screenshot:
% xlim(gca,[0 1]); ylim(gca,[ymin ymax]);         % main axes
% For the insets, un-comment inside createfigure_Flux:
%   xlim(axes2,[0 0.01]); ylim(axes2,[ymin ymax]);
%   xlim(axes3,[0 1]);    ylim(axes3,[-3 3]);
