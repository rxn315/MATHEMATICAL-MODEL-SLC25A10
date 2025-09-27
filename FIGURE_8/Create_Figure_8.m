% ---------------- Call your figure function ----------------
% Inputs: X1, Y1..Y9 in your desired order
% a) [M_m], b) [M_ims], c) Mal. Gradient, d) [S_m], e) [S_ims],
% f) Succ. Gradient, g) [P_m], h) [P_ims], i) Pho. Gradient
run('run_FIG8.m')

createfigure(Time, ...
             Mm, Mims, Mal_grad, ...
             Sm, Sims, Suc_grad, ...
             Pm, Pims, Pho_grad);
