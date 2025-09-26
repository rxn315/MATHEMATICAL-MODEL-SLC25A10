
run('run_Thermo_Validation')
% X axis
X1 = Time(:);

% Build the Y matrices (ensure column vectors)
YMatrix1 = [J_S(:),  J_P(:)];                     % a) fluxes (Jsucc, Jpho)
YMatrix2 = [delta_G(:),  zeros(numel(Time),1)];   % b) ΔG and zero line
YMatrix3 = [Sm(:),   Succ_m_eq*ones(numel(Time),1)];   % c) Sm
YMatrix4 = [Sims(:), Succ_ims_eq*ones(numel(Time),1)]; % d) Sims
YMatrix5 = [Pm(:),   P_m_eq*ones(numel(Time),1)];      % e) Pm
YMatrix6 = [Pims(:), P_ims_eq*ones(numel(Time),1)];    % f) Pims

% Render the figure
createfigure(X1, YMatrix1, YMatrix2, YMatrix3, YMatrix4, YMatrix5, YMatrix6);

