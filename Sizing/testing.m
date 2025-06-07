clear;clc;clf







asdf = struct( ...
    "asdf", 2, ...
    "qwer", 3 ...
    );


station_results = repmat(struct(), 9, 1);

station_results(1) = asdf

























% % Things we have after station 1
% % C_theta_1
% % Updated annulus geometry
% % T01
% % P01
% % T1
% % P1
% % H01
% 
% 
% % Things to calculate for station 1
% % T02 
% 
% 
% %% Station 1 to Station 2 Linking
% % Initialization fluff
% U2 = ang_vel.*r_grid{2};
% U_rotor = (U1 + U2) ./ 2;
% 
% W_m_2 = W_m_1;      % Initial assumption
% C_theta_2 = eleven_five(C_theta_1, psi_c_vec(1), U_r_mean_vec(1), r_mean_vec(1), r_grid{1}, r_grid{2});
% 
% 
% % Transition Thermodynamics
% dH0 = phi_c_vec(1) .* U_rotor.^2;
% T02 = T01 + dH0./cp;
% 

% 
% function [T02, P02, h02] = transition_calcs(H01, U_1, U_2, C_theta_1, W_m_1, cp, R, T01, P01)     %(W_m, C_theta_1, U_1, U_rotor, S1, T01, P01, cp, gamma, loss, phi_c_vec, R)
%     I1 = H01 - U_1.*C_theta_1;
%     W_theta_1 = C_theta_1 - U_1;
%     W_1 = sqrt(W_m_1.^2 + W_theta_1.^2);
%     H01R = I1 + 0.5.*U_1.^2;
% 
% 
%     P01R = "huhhhh"
%     T01R = "huhhhh"
% 
%     I2 = I1;
%     H02R = I2 + 0.5.*U_2.^2;
% 
%     P02R_i = P01R;
%     P02R = P02R_i - loss .* (P01R - P1);
%     T02R = T02R;
% 
%     % C_m = W_m;
%     % C_1 = sqrt(C_m.^2 + C_theta_1.^2);
%     % 
%     % dH0 = phi_c_vec(1) .* U_rotor.^2;
%     % T02 = T01 + dH0./cp;
%     % 
%     % h02 = cp.*T02;                              % Total enthalpy        | Assumed constant spanwise
%     % P2 = P02 .* (T2./T02).^(gamma/(gamma-1));   % Static pressure       | Spanwise distribution
% end





% Turbofan = struct( ...
%     "Info", [],...
%     "Compressor", struct( ...
%         "HP", repmat(struct( ...
%             "T01", [], ...
%             "P01", []), ...
%             num_HP_compressors, 1), ...
%         "LP", repmat(struct( ...
%             "T01", [], ...
%             "P01", []), ...
%             num_LP_compressors, 1) ...
%         ), ...
%     "Combustor", [], ...
%     "Turbine", struct( ...
%         "HP", repmat(struct( ...
%             "T01", [], ...
%             "P01", []), ...
%             num_HP_turbines, 1), ...
%         "LP", repmat(struct( ...
%             "T01", [], ...
%             "P01", []), ...
%             num_LP_turbines, 1) ...
%         ) ...
%     );