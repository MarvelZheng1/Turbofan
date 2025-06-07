clear;clc;clf;
load("Turbofan.mat")

%% Design Parameters
P31_ratio = 1.5;                            % Pressure ratio, P02/P01
eta_c = Turbofan.Specs.Efficiencies.c_lp;   % Stage efficiency
gamma = Turbofan.Specs.Gammas.c_lp;         % Specific-heat ratio
Mrel_max = 0.75;                            % Allowable relative mach number


%% Initial Selections (To be optimized)
J = 1;          % Tip speed parameter, U_tip/a01
zeta = 0.6;     % Hub-to-tip ratio rh/rt
B = -0.4;       % Swirl parameter B

%% Geometry Initialization







x = linspace(1, 1.4, 25);
y = ones(25,3);

for i = 1:25
    y(i,1) = preliminary_design(x(i), eta_c, gamma, Mrel_max, J, zeta, -0.35).m_dot_thing;
    y(i,2) = preliminary_design(x(i), eta_c, gamma, Mrel_max, J, zeta, -0.40).m_dot_thing;
    y(i,3) = preliminary_design(x(i), eta_c, gamma, Mrel_max, J, zeta, -0.45).m_dot_thing;
end

% plot(x,y)

function prelim_stage = preliminary_design(P31_ratio, eta_c, gamma, Mrel_max, J, zeta, B)
    y_vec = linspace(zeta, 1, 50);
    
    tip_work_ratio      = IV_6    (P31_ratio, gamma, eta_c, J);
    rm_swirl_ratio      = IV_18   (zeta, tip_work_ratio);
    [rt_swirl_ratio, A] = IV_S3   (rm_swirl_ratio, zeta, B);
    tip_Cz_ratio        = IV_10   (Mrel_max, J, gamma, rt_swirl_ratio);
    s1_Cz_ratios        = IV_15   (y_vec, A, B, tip_Cz_ratio);
    s2_Cz_ratios        = IV_16   (y_vec, A, B, tip_Cz_ratio, tip_work_ratio, zeta);
    s1_Ctheta_ratios    = IV_S5_1 (y_vec, A, B);
    s2_Ctheta_ratios    = IV_S5_2 (y_vec, A, B, tip_work_ratio);
    s1_rho_ratios       = IV_3    (y_vec, s1_Ctheta_ratios, s1_Cz_ratios, gamma, J);
    m_dot_thing         = IV_2    (y_vec, s1_Cz_ratios, s1_rho_ratios, gamma, J);
    
    prelim_stage = struct( ...
        "params", struct( ...
            "P31_ratio", P31_ratio, ...
            "eta_c",     eta_c, ...
            "gamma",     gamma, ...
            "Mrel_max",  Mrel_max, ...
            "J",         J, ...
            "zeta",      zeta, ...
            "A",         A, ...
            "B",         B ...
            ), ...
        "ratio", struct( ...
            "tip_work_ratio", tip_work_ratio, ...
            "rm_swirl_ratio", rm_swirl_ratio, ...
            "rt_swirl_ratio", rt_swirl_ratio, ...
            "tip_Cz_ratio",   tip_Cz_ratio ...
            ), ...
        "ratios", struct( ...
            "s1_Cz_ratios",     s1_Cz_ratios, ...
            "s2_Cz_ratios",     s2_Cz_ratios, ...
            "s1_Ctheta_ratios", s1_Ctheta_ratios, ...
            "s2_Ctheta_ratios", s2_Ctheta_ratios, ...
            "s1_rho_ratios",    s1_rho_ratios ...
            ), ...
        "m_dot_thing", m_dot_thing ...
        );

    function [rt_swirl_ratio, A] = IV_S3(rm_swirl_ratio, zeta, B)
        ym = (1+zeta)/2;
        A = ym*(rm_swirl_ratio-(B/(ym^2)));
        rt_swirl_ratio = A+B;
    end
    function s1_Ctheta_ratios    = IV_S5_1(y, A, B)
        s1_Ctheta_ratios = A./y + B./(y.^2);
    end
    function s2_Ctheta_ratios    = IV_S5_2(y, A, B, tip_work_ratio)
        s2_Ctheta_ratios = A./y + (B+tip_work_ratio)./(y.^2);
    end
    
    function m_dot_thing    = IV_2 (y, s1_Cz_ratios, s1_rho_ratios, gamma, J)
        constant = pi/2*sqrt(gamma)*J;
        integral = 0;
        for i = 1:length(y)-1
            dy = y(i+1)-y(i);
            integral = integral + s1_rho_ratios(i) * s1_Cz_ratios(i) * y(i) * dy;
        end
        m_dot_thing = constant*integral;
    end
    function s1_rho_ratios  = IV_3 (y, s1_Ctheta_ratios, s1_Cz_ratios, gamma, J)
        s1_rho_ratios = (1 - ((gamma-1)/2*J^2).*(y.^2).*(s1_Ctheta_ratios.^2 + s1_Cz_ratios.^2)).^(gamma/(gamma-1));
    end
    function tip_work_ratio = IV_6 (P31_ratio, gamma, eta_c, J)
        tip_work_ratio = (P31_ratio^((gamma-1)/gamma) - 1) / ((gamma-1)*eta_c*(J.^2));
    end
    function tip_Cz_ratio   = IV_10(Mrel_max, J, gamma, rt_swirl_ratio)
        tip_Cz_ratio = sqrt((Mrel_max^2 * (((1/J)^2) - ((gamma-1)/2) * (rt_swirl_ratio^2)) - (1 - rt_swirl_ratio)^2) / (1 + ((gamma-1)/2) * Mrel_max^2));
    end
    function s1_Cz_ratios   = IV_15(y, A, B, tip_Cz_ratio)
        s1_Cz_ratios = 1./y .* sqrt(tip_Cz_ratio^2 - (2*A^2).*log(y) - (2*A*B).*(1-1./y));
    end
    function s2_Cz_ratios   = IV_16(y, A, B, tip_Cz_ratio, tip_work_ratio, zeta)
        ym = (1+zeta)/2;
        s2_Cz_ratios = 1./y .* sqrt(tip_Cz_ratio^2*ym^2 - (2*A^2).*log(y./ym) - (2*A*(B+tip_work_ratio)).*(1./ym-1./y));
    end
    function rm_swirl_ratio = IV_18(zeta, tip_work_ratio)
        ym = (1+zeta)/2;
        rm_swirl_ratio = 0.5*(1 - tip_work_ratio/(ym^2));
    end
end