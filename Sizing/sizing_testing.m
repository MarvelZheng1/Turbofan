clear;clc;clf

M_f = 0.85;                     % Flight mach number 
% g_cold = 1.401;               % Pre-combustion specific heat ratio (assumed constant)
% g_hot = 1.33;                 % Post-combustion specific heat ratio (assumed constant) 
% specHeat_cold = 1004;         % specific heat for cold components [J/kgK]
% specHeat_hot = 1156;          % specific heat for hot components [J/kgK]
QR = 45000000;                  % heat of reaction of jetA [J/kg]

target_thrust = 150 * 4.44822;   % Newtons (converted from lbf by *4.44822)

combustion_temp = 1233.15;         % Temperature of combustion (burner outlet temperature, T04)    | Kelvin
Rp = 287;                       % Gas constants of products
Ra = 287;                       % Gas constants of air

bypass = 3.3;                     % Bypass ratio!! (fan diameter = 2*(core diameter)

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
    "f",    1.5, ...        % Fan           outlet/inlet, P02/P015
    "c", 10, ...        % Compressor    outlet/inlet, P03/P02
    "b",    1  ...       % Burner         outlet/inlet, P04/P03
    );

% Efficiencies -> 1 = 100%
eta = struct( ...
    "d",     0.94, ...      % Diffuser
    "f",     0.85, ...      % Fan
    "fn",    0.98, ...      % Fan Nozzle
    "c",  0.83, ...      % LP Compressor
    "b",     1.00, ...      % Burner
    "t",  0.89, ...      % LP Turbine
    "n",     0.98  ...      % Nozzle
    );

% Specific Heat Ratios (gamma)
g = struct( ...
    "a",     1.400, ...     % Ambient
    "d",     1.400, ...     % Diffuser
    "f",     1.400, ...     % Fan
    "fn",    1.400, ...     % Fan Nozzle
    "c",  1.400, ...     % LP Compressor
    "b",     1.300, ...     % Burner
    "t",  1.320, ...     % LP Turbine
    "n",     1.340  ...     % Nozzle
    );

%% Ambient Conditions (v=0, so static = total)
Ta = 298;             % Ambient static/total temperature    | Kelvin
Pa = 101300;          % Ambient static/total pressure       | Pascals

u_a = M_f*sqrt(g.a*Ra*Ta);

%% lmao who cares about station 1 am i right (???)
%% Station 1.5: Diffuser Outlet / Fan Inlet / Compressor Inlet
T015 = Ta*(1 + (g.a-1)/2 * M_f^2);
P015 = Pa*(1 + eta.d*(T015/Ta - 1))^(g.d/(g.d-1));

%% Station 2: Fan Outlet
Cpf = (Ra*g.f)/(g.f-1);     % Specific heat of fan

T02 = T015*(1 + 1/eta.f*(Pr.f^((g.f-1)/g.f)-1));
P02 = P015*Pr.f;

%% Station 3: Compressor Outlet/burner Inlet
Cpc = (Ra*g.c)/(g.c-1);     % Specific heat of compressor

T03 = T015*(1 + 1/eta.c*(Pr.c^((g.c-1)/g.c)-1));        % bruh moment
P03 = P015*Pr.c;                                        % bruh moment


%% Station 4: Burner Outlet/Turbine Inlet
Cpb = (Rp*g.b)/(g.b-1);     % Specific heat of burner

T04 = combustion_temp;
P04 = P03*Pr.b;

fr = (T04/T03 - 1) / ((eta.b*QR)/(Cpb*T03)-T04/T03);   % Fuel-air ratio


%% Station 5: Turbine Outlet/Nozzle Inlet
Cpt = (Rp*g.t)/(g.t-1);     % Specific heat of turbine

T05 = ((1+fr)*T04*Cpt - Cpc*(T03-T015) - bypass*Cpf*(T02-T015)) / ((1+fr)*Cpt);
P05 = P04*(1 - 1/eta.t*(1 - T05/T04))^(g.t/(g.t-1));

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

fprintf("ST: %3f\n", ST)
fprintf("TSFC: %3f\n", TSFC)
fprintf("np: %3f\n", eta_p)
fprintf("nth: %3f\n", eta_th)
fprintf("n0: %3f\n", eta_0)
