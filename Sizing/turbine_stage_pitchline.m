% Calculates velocity triangles for a turbine stage at pitchline
function tst = turbine_stage_pitchline(Mc_2m, alpha_1m, alpha_2m, T0_1m, P0_1m, r_mean, ang_vel, gamma_t, R_t, Cp_t, m_dot_t)
    T_2m = T0_1m/(1 + (gamma_t-1)/2*Mc_2m^2);
    a_2m = sqrt(gamma_t*R_t*T_2m);
    C_2m = Mc_2m * a_2m;
    
    Ctheta_2m = C_2m * sind(alpha_2m);
    z_2m = C_2m * cosd(alpha_2m);
    
    z_1m = z_2m;    % assume constant axial velocity
    z_3m = z_2m;    % assume constant axial velocity
    C_1m = z_2m/cosd(alpha_1m);
    Ctheta_1m = C_1m * sind(alpha_1m);
    
    T0_2m = T0_1m;
    T_1m = T0_2m - C_1m^2/(2*Cp_t);
    a_1m = sqrt(gamma_t*R_t*T_1m);
    Mc_1m = C_1m/a_1m;
    Mz_1m = z_1m/a_1m;
    
    zweifel = 1;
    solidityZweifel_s = 2*cosd(alpha_2m)/cosd(alpha_1m)*sind(alpha_2m - alpha_1m);
    
    Ctheta_mean = (Ctheta_1m + Ctheta_2m)/2;  % Initial approximation of average swirl, assuming swirl-free exit
    alpha_2mean = atand(Ctheta_mean/z_2m);
    % optimalSolidity_s = solidityZweifel_s/cosd(alpha_2mean);
    optimalSolidity_s = 2;
    
    o_s_ratio = cosd(alpha_2m)/AAstar_ratio(Mc_2m, gamma_t);
    
    U_2m = ang_vel * r_mean;
    U_1m = U_2m;
    U_3m = U_2m;
    
    Wtheta_2m = Ctheta_2m - U_2m;
    Wtheta_1m = Ctheta_1m - U_2m;
    
    Mw_3m = 0.8;
    Wtheta_3m = -sqrt( (Mw_3m^2*(a_2m^2+(gamma_t-1)*Wtheta_2m^2/2)-z_2m^2) / (1+(gamma_t-1)*Mw_3m^2/2) );
    Ctheta_3m = U_2m + Wtheta_3m;
    W_2m = sqrt(z_2m^2 + Wtheta_2m^2);
    W_3m = sqrt(z_3m^2 + Wtheta_3m^2);
    C_3m = sqrt(z_3m^2 + Ctheta_3m^2);
    W_1m = sqrt(z_1m^2 + Wtheta_1m^2);
    
    beta_1m = acosd(z_1m/W_1m);
    beta_2m = acosd(z_2m/W_2m);
    beta_3m = -acosd(z_3m/W_3m);
    alpha_3m = atand(Ctheta_3m/z_3m);
    
    a_3m = W_3m/Mw_3m;
    Mw_2m = W_2m/a_2m;
    
    degR_1m = 1 - (Ctheta_2m + Ctheta_3m)/(2*U_2m);
    
    T0_2Rm = T_2m + W_2m^2/(2*Cp_t);
    
    profileLoss_s = 0.06;
    
    P_1m = P0_1m * (1 + (gamma_t-1)/2*Mc_1m^2)^(-gamma_t/(gamma_t-1));
    P0_2m = -profileLoss_s*(P0_1m - P_1m)+P0_1m;
    
    P_2m = P0_2m * (1 + (gamma_t-1)/2*Mc_2m^2)^(-gamma_t/(gamma_t-1));
    P0_2Rm = P_2m * (1 + (gamma_t-1)/2*Mw_2m^2)^(gamma_t/(gamma_t-1));
    T0_3m = T0_2m + U_2m*(Ctheta_3m-Ctheta_2m)/Cp_t;
    T_3m = T0_3m - C_3m^2/(2*Cp_t);
    T0_3Rm = T_3m + W_3m^2/(2*Cp_t);
    a_3m = sqrt(gamma_t*R_t*T_3m);
    Mc_3m = C_3m/a_3m;
    
    profileLoss_r = 0.08;
    P0_3Rm = -profileLoss_r*(P0_2Rm - P_2m)+P0_2Rm;
    P_3m = P0_3Rm * (1 + (gamma_t-1)/2*Mw_3m^2)^(-gamma_t/(gamma_t-1));
    P0_3m = P_3m * (1 + (gamma_t-1)/2*Mc_3m^2)^(gamma_t/(gamma_t-1));
    
    % Rotor solidity
    optimalZweifel = 1.0;
    solidityZweifel_r = 2*cosd(beta_3m)/cosd(beta_2m)*sind(beta_2m - beta_3m)/optimalZweifel;
    Wtheta_mean = (Wtheta_2m + Wtheta_3m)/2;
    beta_3mean = atand(Wtheta_mean/z_3m);
    optimalSolidity_r = solidityZweifel_r/cosd(beta_3mean);
    
    % Rotor deviation angle
    devAng = (beta_3m - beta_2m)/(8*optimalSolidity_r);
    
    % Power
    P_spec = U_2m*(Ctheta_2m-Ctheta_3m);
    P = P_spec * m_dot_t;

    % Return "turbine stage triangles"

    tst = struct( ...
        "C_1m", C_1m, ...
        "C_2m", C_2m, ...
        "C_3m", C_3m, ...
        "W_1m", W_1m, ...
        "W_2m", W_2m, ...
        "W_3m", W_3m, ...
        "U_1m", U_1m, ...
        "U_2m", U_2m, ...
        "U_3m", U_3m, ...
        "z_1m", z_1m, ...
        "z_2m", z_2m, ...
        "z_3m", z_3m, ...
        ...
        "Ctheta_1m", Ctheta_1m, ...
        "Ctheta_2m", Ctheta_2m, ...
        "Ctheta_3m", Ctheta_3m, ...
        "Wtheta_1m", Wtheta_1m, ...
        "Wtheta_2m", Wtheta_2m, ...
        "Wtheta_3m", Wtheta_3m, ...
        ...
        "alpha_1m", alpha_1m, ...
        "alpha_2m", alpha_2m, ...
        "alpha_3m", alpha_3m, ...
        "beta_1m", beta_1m, ...
        "beta_2m", beta_2m, ...
        "beta_3m", beta_3m, ...
        ...
        "degR_1m", degR_1m, ...
        "T0_1m", T0_1m, ...
        "T0_2m", T0_2m, ...
        "T0_3m", T0_3m, ...
        "P0_1m", P0_1m, ...
        "P0_2m", P0_2m, ...
        "P0_3m", P0_3m ...
        );
    
    function AAstar = AAstar_ratio(M, gamma)
        AAstar = ((gamma+1)/2) ^ (-(gamma + 1)/(2*(gamma-1))) * (1 + (gamma-1)/2*M^2) ^ ((gamma + 1)/(2*(gamma-1))) / M ;
    end
end