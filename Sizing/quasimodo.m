function [W_m, r_hub_vec] = quasimodo(stage_num, stage_pos, target_m_dot, P0, T0, )
    % ======== Inputs ========
    % stage_num             number of the stage
    % stage_pos             1 = rotor inlet, 2 = stator inlet, 3 = stator outlet
    
    
    % ==== Wm_init Initialization ====
    Wm_init_1 = phi_c_vec(1) * U_r_mean_vec(1);
    
    current_m_dot = target_m_dot + 1;
    new_radius_ish = r_hub_vec(1);
    while abs(current_m_dot - target_m_dot) > 0.001
        r_hub_vec(1:end) = new_radius_ish;
        [r_grid, y_grid, r_mean_vec, U_r_mean_vec, mean_index] = set_grid(r_shroud_vec, r_hub_vec, num_stations, num_surfaces, ang_vel);
    
        C_theta_1 = eleven_four(R_c_vec(1), psi_c_vec(1), U_r_mean_vec(1), r_mean_vec(1), r_grid{1}, 1, 1);
        U_1 = ang_vel .* r_grid{1};
        W_theta_1 = C_theta_1 - U_1;
    
        C_1 = sqrt(Wm_init_1.^2 + C_theta_1.^2);
        h01 = cp*T01;                              % Total enthalpy        | Assumed constant spanwise
    
        L_1 = ones(size(U_1)).*0.03;             % Entropy               | Spanwise distribution        | Constant every rotation (for now)
        I_1 = h01 - U_1.*C_theta_1;              % Rothalpy              | Spanwise distribution        | Changes every rotation  (dependent on C_theta_1
        S_1 = zeros(size(U_1));                  % Entropy               | Spanwise distribution        | Constant every rotation (for now, depends on L_1??)
        T_1 = T01 - C_1.^2./(2*cp);              % Static temperature    | Spanwise distribution        | Changes evert rotation  (dependent on C_!
    
        W_m_1 = seven_fifteen(y_grid{1}, r_grid{1}, C_theta_1, W_theta_1, I_1, S_1, T_1, Wm_init_1, mean_index);
        rho = rho_distro(W_m_1, C_theta_1, T01, P01, cp, gamma, R);
        [new_radius, new_radius_ish, current_m_dot] = annulus_adjust(y_grid{1}, r_grid{1}, rho, W_m_1, target_m_dot);
        fprintf("->")
    end
end