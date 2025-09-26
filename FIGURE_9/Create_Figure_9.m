% ---------------- Call your figure function ----------------
% If you saved it as createfigure_flux.m, change the function name below.
run('run_FIG9.m')

createfigure( ...
    YData_efflux, XData1, ...  % patch 1
    YData_influx, X1, ...      % patch 2 + time for lines
    YMatrix1, X2, Y1);         % 3 fluxes + total

