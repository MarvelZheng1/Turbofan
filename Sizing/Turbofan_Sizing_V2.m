clear;clc;

kelvin_convert = 273;

%% Cruise Conditions
M_0 = 0.88;     % M
P_0 = 15000;    % Pa
T_0 = -40;      % degC

gamma_c = 1.4;
gamma_t = 1.33;

Cp_c = 1004;
Cp_t = 1152;


pi_d = 0.995;   % diffuser
pi_f = 1.6;     % fan
pi_fn = 0.95;   % fan nozzle, convergent
pi_c = 40;      % compressor
pi_b = 0.95;    % burner
pi_n = 0.98;    % main nozzle, convergent

eta_b = 0.992;
eta_m = 0.95;

tau_lambda = 8;
Q_R = 42000000; % J

e_f = 0.9;
e_c = 0.9;
e_t = 0.85;

alpha = 8;

%% Station 0: Flight Conditions
T_0 = T_0 + kelvin_convert;
a_0 = sqrt((gamma_c - 1)*Cp_c*T_0);
V_0 = M_0 * a_0;
P0_0 =  P_0 * (1 + (gamma_c - 1) * M_0^2/2)^(gamma_c/(gamma_c-1));
T0_0 =  T_0 * (1 + (gamma_c - 1) * M_0^2/2);

%% Station 2: Fan Inlet Face
P0_2 = P0_0 * pi_d;
T0_2 = T0_0;

%% Station 13: Fan Outlet
P0_13 = P0_2 * pi_f;
tau_f = pi_f^((gamma_c - 1)/(e_f*gamma_c));
T0_13 = T0_2 * tau_f;

%% Station 19: Fan Nozzle Outlet
P0_19 = P0_13 * pi_fn;
T0_19 = T0_13;
T_19 = T0_19 / (1+(gamma_c-1)/2*)































