clear;clc;

load("Turbofan.mat")

%% Design Variables
m_dot_f = 3;
m_dot_t = 93;     % mass flow, kg/s, air plus fuel
m_dot_c = 100;  % compressor mass flow, just air
rpm = 8770;     % rpm
ang_vel = rpm * 2*pi / 60;

T0_2c = 364.987;
T0_3c = 771;
T0_1m = 1500;    % turbine inlet total temperature, K, is T0_4 in cycle station numbers
P0_1m = 1500000; % turbine inlet total pressure, Pa    is P0_4 in cycle station numbers

r_mean_c = 0.4; % Comrpessor pitchline radius, meters

m_dot_cool = 10;    % Cooling air bleedoff mass flow, kg/s
T0_cool = 704;      % Cooling air temperature, kelvin
P0_cool = 1500000;  % Cooling air pressure, Pa

M_0 = 0;
P0_0m = 101000;
T0_0m = 289.15;

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
T0_5m = double(solve(power_t == m_dot_c*((1-ep)*Cp_t*(T0_1m-T0_5_sym) + ep*Cp_c*(T0_cool-T0_5_sym)), T0_5_sym));
% T0_5m = T0_4 - power_t/(m_dot_t*Cp_t);
deltaT = T0_1m-T0_5m;

Mc_2m = 1.1;
alpha_1m = 0; % assume purely axial inlet velocity
alpha_2m = 60;

%% Here we go fr this time
stage1_tst = turbine_stage_pitchline(Mc_2m, alpha_1m, alpha_2m, T0_1m, P0_1m, r_mean_c, ang_vel, gamma_t, R_t, Cp_t, m_dot_t)

Mc_2m = 0.964;
alpha_2m = 53;
r_mean = 0.42;
stage2_tst = turbine_stage_pitchline(Mc_2m, stage1_tst.alpha_3m, alpha_2m, stage1_tst.T0_3m, stage1_tst.P0_3m, r_mean, ang_vel, gamma_t, R_t, Cp_t, m_dot_t)