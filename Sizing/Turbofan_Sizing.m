clear;clc;clf

M_f = 0.85;                     % Flight mach number 
target_thrust = 5000 * 4.44822;   % Newtons (converted from lbf by *4.44822)
bypass = 3;                     % Bypass ratio!! (3: fan diameter = 2*(core diameter)

QR = 45000000;                  % heat of reaction of jetA [J/kg]
combustion_temp = 1750;         % Temperature of combustion (burner outlet temperature, T04)    | Kelvin
Rp = 287;                       % Gas constants of products
Ra = 287;                       % Gas constants of air


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
    "cLP", sqrt(5), ...        % LP Compressor    outlet/inlet, P03/P02
    "cHP", sqrt(5), ...        % HP Compressor    outlet/inlet, P03/P02
    "b",    1  ...       % Burner         outlet/inlet, P04/P03
    );

% Efficiencies -> 1 = 100%
eta = struct( ...
    "d",     0.94, ...      % Diffuser
    "f",     0.85, ...      % Fan
    "fn",    0.98, ...      % Fan Nozzle
    "cLP",  0.83, ...      % LP Compressor
    "cHP",  0.83, ...      % HP Compressor
    "b",     1.00, ...      % Burner
    "tHP",  0.89, ...      % HP Turbine
    "tLP",  0.89, ...      % LP Turbine
    "n",     0.98  ...      % Nozzle
    );

% Specific Heat Ratios (gamma)
g = struct( ...
    "a",     1.400, ...     % Ambient
    "d",     1.400, ...     % Diffuser
    "f",     1.400, ...     % Fan
    "fn",    1.400, ...     % Fan Nozzle
    "cLP",  1.400, ...     % LP Compressor
    "cHP",  1.400, ...     % HP Compressor
    "b",     1.300, ...     % Burner
    "tHP",  1.320, ...     % HP Turbine
    "tLP",  1.320, ...     % LP Turbine
    "n",     1.340  ...     % Nozzle
    );

%% Ambient Conditions (v=0, so static = total)
T_0 = 298;             % Ambient static/total temperature    | Kelvin
P_0 = 101300;          % Ambient static/total pressure       | Pascals

u_a = M_f*sqrt(g.a*Ra*T_0);

%% lmao who cares about station 1 am i right (???)
%% Station 1.5: Diffuser Outlet /Fan Inlet
T0_15 = T_0*(1 + (g.a-1)/2 * M_f^2);
P0_15 = P_0*(1 + eta.d*(T0_15/T_0 - 1))^(g.d/(g.d-1));

%% Station 2: Fan Outlet/LP Compressor Inlet
Cp_f = (Ra*g.f)/(g.f-1);     % Specific heat of fan

T0_2 = T0_15*(1 + 1/eta.f*(Pr.f^((g.f-1)/g.f)-1));
P0_2 = P0_15*Pr.f;

%% Station 2.5: LP Compressor Outlet/HP Compressor Inlet
Cp_cLP = (Ra*g.cLP)/(g.cLP-1);     % Specific heat of compressor

T0_25 = T0_2*(1 + 1/eta.cLP*(Pr.cLP^((g.cLP-1)/g.cLP)-1));
P0_25 = P0_2*Pr.cLP;

%% Station 3: HP Compressor Outlet/Burner Inlet
Cp_cHP = (Ra*g.cHP)/(g.cHP-1);     % Specific heat of compressor

T0_3 = T0_25*(1 + 1/eta.cHP*(Pr.cHP^((g.cHP-1)/g.cHP)-1));
P0_3 = P0_25*Pr.cHP;

%% Station 4: Burner Outlet/HP Turbine Inlet
Cp_b = (Rp*g.b)/(g.b-1);     % Specific heat of burner

T0_4 = combustion_temp;
P0_4 = P0_3*Pr.b;

fr = (T0_4/T0_3 - 1)/((eta.b*QR)/(Cp_b*T0_3)-T0_4/T0_3);   % Fuel-air ratio

%% Station 4.5: HP Turbine Outlet/LP Turbine Inlet
Cp_tHP = (Rp*g.tHP)/(g.tHP-1);     % Specific heat of turbine

T0_45 = ((1+fr)*T0_4*Cp_tHP - Cp_cHP*(T0_3-T0_25)) / ((1+fr)*Cp_tHP);
P0_45 = P0_4*(1 - 1/eta.tHP*(1 - T0_45/T0_4))^(g.tHP/(g.tHP-1));

%% Station 5: LP Turbine Outlet/Nozzle Inlet
Cp_tLP = (Rp*g.tLP)/(g.tLP-1);     % Specific heat of turbine

T0_5 = ((1+fr)*T0_45*Cp_tLP - Cp_cLP*(T0_25-T0_2) - bypass*Cp_f*(T0_2-T0_15)) / ((1+fr)*Cp_tLP);
P0_5 = P0_45*(1 - 1/eta.tLP*(1 - T0_5/T0_45))^(g.tLP/(g.tLP-1));

%% Station 6: Afterburner (there is none lmao)
T0_6 = T0_5;
P0_6 = P0_5;

%% Station 7: Nozzle Outlet
T0_7 = T0_6;
P0_7 = P0_6;

%% Station 8: Aft Ambient
T_8 = T_0;
P_8 = P_0;

%% Nozzle Exit Velocities
u_ec = sqrt(2*eta.n *(g.n /(g.n -1))*Rp*T0_7*(1 - (P_8/P0_7)^((g.n -1)/g.n )));
u_ef = sqrt(2*eta.fn*(g.fn/(g.fn-1))*Ra*T0_2*(1 - (P_8/P0_2)^((g.fn-1)/g.fn)));

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
Info.core_mass_flow_air = Info.mass_flow_air/(1+bypass);

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
        "T0", T_0, ...
        "P0", P_0), ...
    "S15", struct( ...
        "T0", T0_15, ...
        "P0", P0_15), ...
    "S2", struct( ...
        "T0", T0_2, ...
        "P0", P0_2), ...
    "S25", struct( ...
        "T0", T0_25, ...
        "P0", P0_25), ...
    "S3", struct( ...
        "T0", T0_3, ...
        "P0", P0_3), ...
    "S4", struct( ...
        "T0", T0_4, ...
        "P0", P0_4), ...
    "S45", struct( ...
        "T0", T0_45, ...
        "P0", P0_45), ...
    "S5", struct( ...
        "T0", T0_5, ...
        "P0", P0_5), ...
    "S6", struct( ...
        "T0", T0_6, ...
        "P0", P0_6), ...
    "S7", struct( ...
        "T", T0_7, ...
        "P", P0_7) ...
    );

%% Create Cp Struct
Cp = struct( ...
    "f",    Cp_f, ...
    "cLP", Cp_cLP, ...
    "cHP", Cp_cHP, ...
    "b",    Cp_b, ...
    "tHP", Cp_tHP, ...
    "tLP", Cp_tLP ...
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
plot([0,1.5,2,2.5,3,4,4.5,5,6,7,8], [P_0, P0_15, P0_2, P0_25, P0_3, P0_4, P0_45, P0_5, P0_6, P0_7, P_8])
ylabel("Pressure")
yyaxis right
plot([0,1.5,2,2.5,3,4,4.5,5,6,7,8], [T_0, T0_15, T0_2, T0_25, T0_3, T0_4, T0_45, T0_5, T0_6, T0_7, T_8])
ylabel("Temperature")
xlabel("Station Number")
grid on