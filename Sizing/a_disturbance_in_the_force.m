clear;clc;

load("Turbofan.mat")

%% Design Variables
m_dot_f = 3;
m_dot_t = 93;     % mass flow, kg/s, air plus fuel
m_dot_c = 100;  % compressor mass flow, just air
rpm = 8770;     % rpm

T0_2c = 364.987;
T0_3c = 771;
T0_1 = 1500;    % turbine inlet total temperature, K, is T0_4 in cycle station numbers
P0_1 = 1500000; % turbine inlet total pressure, Pa    is P0_4 in cycle station numbers

r_mean_c = 0.4; % Comrpessor pitchline radius, meters

m_dot_cool = 10;    % Cooling air bleedoff mass flow, kg/s
T0_cool = 704;      % Cooling air temperature, kelvin
P0_cool = 1500000;  % Cooling air pressure, Pa

M_0 = 0;
P0_0 = 101000;
T0_0 = 289.15;

Cp_c = 1004.5;
Cp_t = 1243.67;
gamma_c = 1.40;
gamma_t = 1.30;
R_t = (gamma_t-1)*Cp_t/gamma_t;

eta_m = 0.996;
ep = 0.1;

%% Here we go
power_c = m_dot_c*Cp_c*(T0_3c-T0_2c);
power_t = power_c/eta_m;

syms T0_5_sym;
T0_5m = double(solve(power_t == m_dot_c*((1-ep)*Cp_t*(T0_1-T0_5_sym) + ep*Cp_c*(T0_cool-T0_5_sym)), T0_5_sym));
% T0_5m = T0_4 - power_t/(m_dot_t*Cp_t);
deltaT = T0_1-T0_5m;

Mc_2m = 1.1;
alpha_1m = 0; % assume purely axial inlet velocity
alpha_2m = 60;

T_2m = T0_1/(1 + (gamma_t-1)/2*Mc_2m^2);
a_2m = sqrt(gamma_t*R_t*T_2m);
C_2m = Mc_2m * a_2m;

Ctheta_2m = C_2m * sind(alpha_2m);
z_2m = C_2m * cosd(alpha_2m);

z_1m = z_2m;    % assume constant axial velocity
C_1m = z_2m/cosd(alpha_1m);

T_1m = T0_1 - C_1m^2/(2*Cp_t);
a_1m = sqrt(gamma_t*R_t*T_1m);
Mc_1m = C_1m/a_1m;
Mz_1m = z_1m/a_1m;

zweifel = 1;
solidity_zweifel = 2*cosd(alpha_2m)/cosd(alpha_1m)*sind(alpha_2m - alpha_1m);

Ctheta_mean = Ctheta_2m/2; % Initial approximation of average swirl, assuming swirl-free exit
alpha_2mean = atand(Ctheta_mean/z_2m);
% solidity_optimal = solidity_zweifel/cosd(alpha_2mean);
solidity_optimal = 2;

o_s_ratio = cosd(alpha_2m)/AAstar_ratio(Mc_2m, gamma_t);






function AAstar = AAstar_ratio(M, gamma)
    AAstar = ((gamma+1)/2) ^ (-(gamma + 1)/(2*(gamma-1))) * (1 + (gamma-1)/2*M^2) ^ ((gamma + 1)/(2*(gamma-1))) / M ;
end