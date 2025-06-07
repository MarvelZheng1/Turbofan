%% Incredibly Named Functions
classdef txtbk
properties
end
    methods(Static)
        function C_theta_1 = eleven_four(Rc_c1, psi_c1, U_c1, r_c1, r_1, n, m)
            C_theta_1 = U_c1.*((1-Rc_c1).*(r_c1./r_1).^n - psi_c1./2.*(r_c1./r_1).^m);
        end
        function C_theta_2 = eleven_five(C_theta_1, psi_c1, U_c1, r_c1, r_1, r_2)
            C_theta_2 = (r_1.*C_theta_1 + psi_c1.*U_c1.*r_c1) ./ r_2;
        end
        function C_theta_3 = eleven_six(Rc_c2, psi_c2, U_c3, r_c3, r_3, n, m)
            C_theta_3 = U_c3.*((1-Rc_c2).*(r_c3./r_3).^n - psi_c2./2.*(r_c3./r_3).^m);
        end
    end
end