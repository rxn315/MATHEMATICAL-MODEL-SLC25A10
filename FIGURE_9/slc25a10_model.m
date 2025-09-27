    function dYdt = slc25a10_model(t, y, Vims, Vm, Ks_m, Km_m, Kp_m, Ts_max, Tm_max,lambda21,lambda31)

    Mm = y(1); 
    Mims = y(2); 
    Sm = y(3); 
    Sims = y(4); 
    Pm = y(5); 
    Pims = y(6); 
    

    % -----------------------------
    % Delta Terms
    % -----------------------------
    delta1 = 1 + Sm/Ks_m + Mm/Km_m + Pm/Kp_m;
    delta2 = 1 + Sims/Ks_m + Mims/Km_m + Pims/Kp_m;

    % -----------------------------
    % Phi Terms
    % -----------------------------
    % ----- Phi terms (using lambda ratios) -----
    % phi1: no volume scaling
    phi1_num = Sm + lambda21*Mm + lambda31*Pm;
    phi1_den = Sims + lambda21*Mims + lambda31*Pims;
    phi1 = phi1_num / phi1_den;

    % phi2: phosphate term scaled by (Vm/Vims)^2
    vol2 = (Vm / Vims)^2;
    phi2_num = Sm + lambda21*Mm + vol2*lambda31*Pm;
    phi2_den = Sims + lambda21*Mims + vol2*lambda31*Pims;
    phi2 = phi2_num / phi2_den;
 
    % -----------------------------
    % Flux Equations
    % -----------------------------
    % J_succ
    J_succ = (1/Vm) * (Ts_max * (phi2*Sims - Sm) / (Ks_m * delta1 + Ks_m * phi2 * delta2)) ...
           - (1/Vims) * (Ts_max * (Sm - phi1*Sims) / (Ks_m * delta1 + Ks_m * phi1 * delta2));

    % J_mal
    J_mal = (1/Vm) * (Tm_max * (phi2*Mims - Mm) / (Km_m * delta1 + Km_m * phi2 * delta2)) ...
           - (1/Vims) * (Tm_max * (Mm - phi1*Mims) / (Km_m * delta1 + Km_m * phi1 * delta2));

    % J_pho
    % J_pho = (1/Vims) * (Tp_max * (Pm - phi2*Pims) / (Kp_m * delta1 + Kp_m * phi2 * delta2)) ...
    %        - (1/Vm) * (Tp_max * (phi1*Pims - Pm) / (Kp_m * delta1 + Kp_m * phi1 * delta2));

    % J_pho=(1-exp(-Pm))*(J_succ+J_mal);
    % if Pm>=0 && Pims>=0
    J_pho=(J_succ+J_mal);
    % else
        % J_pho=0;
    % end

   
    % J_pho=J_succ+J_mal;
     % %% Net Cycle
    % % ODEs
    dM_m_dt = J_mal;
    dM_c_dt = -J_mal;
    dS_m_dt = J_succ;
    dS_c_dt = -J_succ;
    dP_m_dt = -J_pho;
    dP_c_dt = J_pho;

    

    % Pack derivatives
    dYdt = [dM_m_dt; dM_c_dt; dS_m_dt; dS_c_dt; dP_m_dt; dP_c_dt];
end


