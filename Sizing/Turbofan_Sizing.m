clear;clc;clf

M_f = 0.85;                     % Flight mach number 
% g_cold = 1.401;               % Pre-combustion specific heat ratio (assumed constant)
% g_hot = 1.33;                 % Post-combustion specific heat ratio (assumed constant) 
% specHeat_cold = 1004;         % specific heat for cold components [J/kgK]
% specHeat_hot = 1156;          % specific heat for hot components [J/kgK]
QR = 45000000;                  % heat of reaction of jetA [J/kg]

target_thrust = 250 * 4.44822;   % Newtons (converted from lbf by *4.44822)

combustion_temp = 1750;         % Temperature of combustion (burner outlet temperature, T04)    | Kelvin
Rp = 287;                       % Gas constants of products
Ra = 287;                       % Gas constants of air

bypass = 3;                     % Bypass ratio!! (fan diameter = 2*(core diameter)

Info = struct( ...
    "M_f",              M_f, ...
    "Rp",               Rp, ...
    "Ra",               Ra, ...
    "combustion_temp",  combustion_temp, ...
    "QR",               QR, ...
    "bypass",           bypass, ...
    "target_thrust",    target_thrust ...
    );

% Design Pressure Ratios
Pr = struct( ...
    "f",    1.2, ...        % Fan           outlet/inlet, P02/P015
    "c_lp", sqrt(4), ...        % Compressor    outlet/inlet, P03/P02
    "c_hp", sqrt(4), ...        % Compressor    outlet/inlet, P03/P02
    "b",    1  ...       % Burner         outlet/inlet, P04/P03
    );

% Efficiencies -> 1 = 100%
eta = struct( ...
    "d",     0.94, ...      % Diffuser
    "f",     0.85, ...      % Fan
    "fn",    0.98, ...      % Fan Nozzle
    "c_lp",  0.83, ...      % LP Compressor
    "c_hp",  0.83, ...      % HP Compressor
    "b",     1.00, ...      % Burner
    "t_hp",  0.89, ...      % HP Turbine
    "t_lp",  0.89, ...      % LP Turbine
    "n",     0.98  ...      % Nozzle
    );

% Specific Heat Ratios (gamma)
g = struct( ...
    "a",     1.400, ...     % Ambient
    "d",     1.400, ...     % Diffuser
    "f",     1.400, ...     % Fan
    "fn",    1.400, ...     % Fan Nozzle
    "c_lp",  1.400, ...     % LP Compressor
    "c_hp",  1.400, ...     % HP Compressor
    "b",     1.300, ...     % Burner
    "t_hp",  1.320, ...     % HP Turbine
    "t_lp",  1.320, ...     % LP Turbine
    "n",     1.340  ...     % Nozzle
    );

%% Ambient Conditions (v=0, so static = total)
Ta = 298;             % Ambient static/total temperature    | Kelvin
Pa = 101300;          % Ambient static/total pressure       | Pascals

u_a = M_f*sqrt(g.a*Ra*Ta);

%% lmao who cares about station 1 am i right (???)
%% Station 1.5: Diffuser Outlet /Fan Inlet
T015 = Ta*(1 + (g.a-1)/2 * M_f^2);
P015 = Pa*(1 + eta.d*(T015/Ta - 1))^(g.d/(g.d-1));

%% Station 2: Fan Outlet/LP Compressor Inlet
Cpf = (Ra*g.f)/(g.f-1);     % Specific heat of fan

T02 = T015*(1 + 1/eta.f*(Pr.f^((g.f-1)/g.f)-1));
P02 = P015*Pr.f;

%% Station 2.5: LP Compressor Outlet/HP Compressor Inlet
Cpc_LP = (Ra*g.c_lp)/(g.c_lp-1);     % Specific heat of compressor

T025 = T02*(1 + 1/eta.c_lp*(Pr.c_lp^((g.c_lp-1)/g.c_lp)-1));
P025 = P02*Pr.c_lp;

%% Station 3: HP Compressor Outlet/Burner Inlet
Cpc_HP = (Ra*g.c_hp)/(g.c_hp-1);     % Specific heat of compressor

T03 = T025*(1 + 1/eta.c_hp*(Pr.c_hp^((g.c_hp-1)/g.c_hp)-1));
P03 = P025*Pr.c_hp;

%% Station 4: Burner Outlet/HP Turbine Inlet
Cpb = (Rp*g.b)/(g.b-1);     % Specific heat of burner

T04 = combustion_temp;
P04 = P03*Pr.b;

fr = (T04/T03 - 1)/((eta.b*QR)/(Cpb*T03)-T04/T03);   % Fuel-air ratio

%% Station 4.5: HP Turbine Outlet/LP Turbine Inlet
Cpt_HP = (Rp*g.t_hp)/(g.t_hp-1);     % Specific heat of turbine

T045 = ((1+fr)*T04*Cpt_HP - Cpc_HP*(T03-T025)) / ((1+fr)*Cpt_HP);
P045 = P04*(1 - 1/eta.t_hp*(1 - T045/T04))^(g.t_hp/(g.t_hp-1));

%% Station 5: LP Turbine Outlet/Nozzle Inlet
Cpt_LP = (Rp*g.t_lp)/(g.t_lp-1);     % Specific heat of turbine

T05 = ((1+fr)*T045*Cpt_LP - Cpc_LP*(T025-T02) - bypass*Cpf*(T02-T015)) / ((1+fr)*Cpt_LP);
P05 = P045*(1 - 1/eta.t_lp*(1 - T05/T045))^(g.t_lp/(g.t_lp-1));

%% Station 6: Afterburner (there is none lmao)
T06 = T05;
P06 = P05;

%% Station 7: Nozzle Outlet/Exit Ambient
T7 = Ta;
P7 = Pa;

%% Nozzle Exit Velocities
u_ec = sqrt(2*eta.n *(g.n /(g.n -1))*Rp*T06*(1 - (Pa/P06)^((g.n -1)/g.n )));
u_ef = sqrt(2*eta.fn*(g.fn/(g.fn-1))*Ra*T02*(1 - (Pa/P02)^((g.fn-1)/g.fn)));

%% Performance Metrics
ST = (1+fr)*u_ec + bypass*u_ef - (1+bypass)*u_a;     % Specific Thrust
TSFC = fr/ST;                                        % Thrust Specific Fuel Consumption
eta_p  = ST*u_a / ((1+fr)*((u_ec^2)/2) + bypass*((u_ef^2)/2) - (1+bypass)*((u_a^2)/2));     % Propulsive Efficiency
eta_th = ((1+fr)*((u_ec^2)/2) + bypass*((u_ef^2)/2) - (1+bypass)*((u_a^2)/2)) / (fr*QR);    % Thermal Efficiency
eta_0  = eta_p*eta_th;                                                                      % Overall Efficiency

Perf = struct( ...
    "ST",           ST, ...
    "TSFC",         TSFC, ...
    "prop_eff",     eta_p, ...
    "thermo_eff",   eta_th, ...
    "overall_eff",  eta_0 ...
    );

%% Mass Flow of Air     | Units: kg/s
syms m_dot
Info.mass_flow = double(vpa(solve(target_thrust == ST * m_dot, m_dot)));
Info.mass_flow_air = Info.mass_flow / (1 + fr);
Info.mass_flow_fuel = Info.mass_flow - Info.mass_flow_air;
Info.core_mass_flow_air = Info.mass_flow_air/bypass;

Info.u_ec = u_ec;
Info.u_ef = u_ef;

fprintf("ST: %3f\n", ST)
fprintf("TSFC: %3f\n", TSFC)
fprintf("np: %3f\n", eta_p)
fprintf("nth: %3f\n", eta_th)
fprintf("n0: %3f\n", eta_0)




%% Create Temperature and Pressure Struct
T0P0 = struct( ...
    "Sa", struct( ...
        "T0", Ta, ...
        "P0", Pa), ...
    "S15", struct( ...
        "T0", T015, ...
        "P0", P015), ...
    "S2", struct( ...
        "T0", T02, ...
        "P0", P02), ...
    "S25", struct( ...
        "T0", T025, ...
        "P0", P025), ...
    "S3", struct( ...
        "T0", T03, ...
        "P0", P03), ...
    "S4", struct( ...
        "T0", T04, ...
        "P0", P04), ...
    "S45", struct( ...
        "T0", T045, ...
        "P0", P045), ...
    "S5", struct( ...
        "T0", T05, ...
        "P0", P05), ...
    "S6", struct( ...
        "T0", T06, ...
        "P0", P06), ...
    "S7", struct( ...
        "T", T7, ...
        "P", P7) ...
    );

%% Create Cp Struct
Cp = struct( ...
    "f",    Cpf, ...
    "c_LP", Cpc_LP, ...
    "c_HP", Cpc_HP, ...
    "b",    Cpb, ...
    "t_HP", Cpt_HP, ...
    "t_LP", Cpt_LP ...
    );

% clear Ta Pa T015 P015 T02 P02 T025 P025 T03 P03 T04 P04 T045 P045 T05 P05 T06 P06 T07 P07

Turbofan = struct( ...
    "Specs", struct( ...
        "Info", Info, ...
        "Efficiencies", eta, ...
        "Gammas", g, ...
        "desPR", Pr ...
        ), ...
    "Thermos", T0P0, ...
    "Cp", Cp, ...
    "Performance", Perf ...
    );

save("Turbofan.mat", "Turbofan")


figure(1)
tiledlayout(1,1, TileSpacing='tight', Padding='tight')
nexttile
title("Temperature and Pressure at Various Stations")
yyaxis left
plot([0,1.5,2,2.5,3,4,4.5,5,6,7], [Pa, P015, P02, P025, P03, P04, P045, P05, P06, P7])
ylabel("Pressure")
yyaxis right
plot([0,1.5,2,2.5,3,4,4.5,5,6,7], [Ta, T015, T02, T025, T03, T04, T045, T05, T06, T7])
ylabel("Temperature")
xlabel("Station Number")
grid on