%CREATEFIGURE_SDH_SCAS_SLC25A10(X1, YMatrix1)
% X1:        time vector (x-axis)
% YMatrix1:  N×3 matrix [J_tot, J_SDH, J_SCAS]
run('Trial.m')
%Example call:
  X1       = Time(:);
  YMatrix1 = [J_tot(:), SDH(Sm,Mm), SCAS(Sm,Pm)];
  createfigure_SDH_SCAS_SLC25A10(Time, YMatrix1);