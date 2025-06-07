clear;clc;close all
%% ======== User-Defined Design Variables ========

rpm = 40000;

num_stages = 4;
num_stations = (num_stages+1) * 2;

num_surfaces = 50;

%% ======== Inputs from the sizing script ========
load("Turbofan.mat")
T01     = ones(1, num_surfaces) * Turbofan.Thermos.S2.T0;
T0_out  = ones(1, num_surfaces) * Turbofan.Thermos.S25.T0;
P01     = ones(1, num_surfaces) * Turbofan.Thermos.S2.P0;
gamma   = Turbofan.Specs.Gammas.c_lp;
R       = Turbofan.Specs.Info.Ra;
cp      = Turbofan.Cp.c_LP;

target_m_dot = Turbofan.Specs.Info.mass_flow_air;

% ======== Dimensionless Performance Parameters, per stage ========
phi_c_vec = [0.6, 0.6, 0.6, 0.6, 0.6];        % Flow coefficient
psi_c_vec = [0.3, 0.3, 0.3, 0.3, 0.3];        % Work coefficient ----> to be determined from temperature rises later
R_c_vec   = [0.5, 0.5, 0.5, 0.5, 0.5];        % Degree of reaction

%% ======== Startup Fluff ========
% Hub and Shroud Geometry Initialization
% Vector along axis
r_shroud_vec = ones(1, num_stations)*0.001* 110;
r_hub_vec    = ones(1, num_stations)*0.001* 75;
% Misc
ang_vel = rpm * 2 * pi / 60;

[r_grid, y_grid, r_mean_vec, U_r_mean_vec, mean_index] = set_grid(r_shroud_vec, r_hub_vec, num_stations, num_surfaces, ang_vel);

[station_results, station_feeder] = create_framework();

station_feeder(1) = struct( ...
    "stage_num",    1, ...      % Stage number
    "stage_123",    1, ...      % Position within stage (local station number)
    "phi_c_vec",    phi_c_vec, ...      % Phi's                                 | Vector per stage
    "psi_c_vec",    psi_c_vec, ...      % Psi's                                 | Vector per stage
    "R_c_vec",      R_c_vec, ...        % Mean radii degrees of reaction        | Vector per absolute station
    "r_hub_vec",    r_hub_vec, ...      % Hub radii                             | Vector per absolute station
    "r_shroud_vec", r_shroud_vec, ...   % Shroud radii                          | Vector per absolute station
    "U_r_mean_vec", U_r_mean_vec, ...   % U at mean radius                      | Vector per absolute station
    "r_mean_vec",   r_mean_vec, ...     % Mean radii                            | Vector per absolute station
    "mean_index",   mean_index, ...     % Index of mean stream surface
    "num_stations", num_stations, ...   % Number of absolute stations
    "num_surfaces", num_surfaces, ...   % Number of stream surfaces
    "ang_vel",      ang_vel, ...`       % Angular velocity
    "T0",           T01, ...             % Total temperature                     | Spanwise vector
    "P0",           P01, ...             % Total pressure                        | Spanwise vector
    "cp",           cp, ...             % Cp
    "gamma",        gamma, ...          % Gamma
    "R",            R, ...              % Gas constant
    "target_m_dot", target_m_dot ...    % Target mass flux
);


%% ======== Platform Nine and Seven Eighths ========
for i = 1:length(station_results)
    fprintf("Absolute station #: %i", i)
    [station_results(i), r_grid, y_grid, station_feeder(i+1)] = stage(station_feeder(i), r_grid, y_grid);
end



fprintf("\n")
for i = 1:length(station_results)
    fprintf("Total Temperature at Station %i: %f\n", i, station_results(i).thermo.T0_cur(station_results(i).mean_index))
end
fprintf("==================================\n")
for i = 1:length(station_results)
    fprintf("Total Pressure at Station %i: %f\n", i, station_results(i).thermo.P0_cur(station_results(i).mean_index))
end
    fprintf("\n==================================\n\n")
for i = 1:length(station_results)
    fprintf("Static Temperature at Station %i: %f\n", i, station_results(i).thermo.T_cur(station_results(i).mean_index))
end
    fprintf("==================================\n")
for i = 1:length(station_results)
    fprintf("Static Pressure at Station %i: %f\n", i, station_results(i).thermo.P_cur(station_results(i).mean_index))
end






% plot_vel_triangle([300,600], triangle1_mean, 1, '-k', '-b', '-r')
% plot_vel_triangle([300,600], triangle2_mean, 1, '-k', '-b', '-r')
% plot_vel_triangle([300,600], triangle3_mean, 1, '-k', '-b', '-r')
% plot_vel_triangle([300,600], triangle4_mean, 1, '-k', '-b', '-r')

figure(Name="Streamlines")
hold on
for i = 1:num_surfaces
    streamline = ones(1,num_stations);
    for j = 1: num_stations
        streamline(j) = r_grid{j}(i);
    end
    plot(streamline, 'k--')
end
plot([station_results(end).r_shroud_vec', station_results(end).r_hub_vec', station_results(end).r_mean_vec'], '.-')
height = max(station_results(end).r_shroud_vec)-min(station_results(end).r_hub_vec);
ylim([min(station_results(end).r_hub_vec)-height/10, max(station_results(end).r_shroud_vec)+height/10])
xlim([0, num_stations+1])

% plot_halftree(000, station1.W_m, 12, -100, 1, '-b')
% plot_halftree(300, station2.W_m, 12, -100, 1, '-b')
% plot_halftree(600, station3.W_m, 12, -100, 1, '-b')
% 
% plot_halftree(000, station1.C_theta, 12, -96, 1, '-r')
% plot_halftree(300, station2.C_theta, 12, -96, 1, '-r')
% plot_halftree(600, station3.C_theta, 12, -96, 1, '-r')




%% Functions
function [T02, P02, T2, P2, T1, P1] = thermobridge(C_theta_1, W_m_1, U1, C_theta_2, W_m_2, U2, T01, P01, cp, rho, loss)

    % ======== C and W Determination ========
    % U1 = U2 = 0 if stator
    C_m_1 = W_m_1;
    C1 = sqrt(C_m_1.^2 + C_theta_1.^2);
    W_theta_1 = C_theta_1 - U1;
    W1 = sqrt(W_m_1.^2 + W_theta_1.^2);

    C_m_2 = W_m_2;
    C2 = sqrt(C_m_2.^2 + C_theta_2.^2);
    W_theta_2 = C_theta_2 - U2;
    W2 = sqrt(W_m_2.^2 + W_theta_2.^2);

    % ======== Station 1 ========
    T1 = T01 - C1.^2./(2.*cp);
    P1 = P01 - rho .* C1.^2 ./ 2;

    T1R = T1;
    P1R = P1;

    T01R = T1R + W1.^2./(2.*cp);
    P01R = P1R + rho.*W1.^2./2;

    % ======== Station 2 ========
    T02Ri = T01R;                       % Assume adiabatic
    P02Ri = P01R;                       % Assume isentropic

    T02R = T02Ri;                       % Assumption becomes reality (O_o)
    P02R = P02Ri - loss.*(P01R - P1R);  % Factor in loss

    T2R = T02R - W2.^2./(2.*cp);
    P2R = P02R - rho .* W2.^2 ./ 2;

    T2 = T2R;
    P2 = P2R;

    T02 = T2 + C2.^2./(2.*cp);
    P02 = P2 + rho.*C2.^2./2;

end




function [r_grid, y_grid, r_mean_vec, U_r_mean_vec, mean_index] = set_grid(r_shroud_vec, r_hub_vec, num_stations, num_surfaces, ang_vel)
    % ==== Mean Radius Initialization ====
    r_mean_vec   = ones(1, num_stations);
    r_mean_real_vec   = ones(1, num_stations);
    
    % ==== y-grid and r_grid Initialization ====
    % y = r - r_hub
    r_grid = cell(1,num_stations); % {1} = leftmost quasi-normal  |  (1) = hub radius  |  Basically {x}(y) is cartesian
    y_grid = cell(1,num_stations); % {1} = leftmost quasi-normal  |  (1) = hub radius  |  Basically {x}(y) is cartesian
    
    % ==== Find meanline radii ====
    syms r_mean_sym;
    for i = 1:num_stations
        r_mean_sols = solve(r_shroud_vec(i)^2-r_mean_sym^2 == r_mean_sym^2-r_hub_vec(i)^2, r_mean_sym);
        r_mean_real_vec(i) = double(r_mean_sols(2));
        r_grid{i} = linspace(r_hub_vec(i), r_shroud_vec(i), num_surfaces);         % Grid of radii values
    end
    
    for i = 1:num_stations
        y_grid{i} = r_grid{i} - r_grid{i}(1);           % Creates y-grid, for meridionial analysis
    end
    
    r_mean_dist = r_grid{1}-r_mean_real_vec(1);
    [~, mean_index] = min(abs(r_mean_dist));
    for i = 1:num_stations
        r_mean_vec(i) = r_grid{i}(mean_index);          % Y-value streamline closest to the mean radius of station 1 (equal annulus area above and below mean radius)
    end

    % Find meanline U velocities (returns mean radius values of U axially)
    U_r_mean_vec = ang_vel.*r_mean_vec;
end

function [station_results, station_feeder] = create_framework()

    results_template = struct( ...
        "r_hub_vec",    [], ...
        "r_shroud_vec", [], ...
        "r_mean_vec",   [], ...
        "U_r_mean_vec", [], ...
        "mean_index",   [], ...
        "W_m",          [], ...
        "C_theta",      [], ...
        "W_theta",      [], ...
        "U",            [], ...
        "Loss",         [], ...
        "rho",          [], ...
        "thermo", struct( ...
            "T0_next",  [], ...
            "P0_next",  [], ...
            "T0_cur",   [], ...
            "P0_cur",   [], ...
            "T_next",   [], ...
            "P_next",   [], ...
            "T_cur",    [], ...
            "P_cur",    [] ...
        ) ...
    );

    station_results = repmat(results_template, 9, 1);

    feeder_template = struct( ...
        "stage_num",    [], ... % Stage number
        "stage_123",    [], ... % Position within stage (local station number)
        "phi_c_vec",    [], ... % Phi's                                 | Vector per stage
        "psi_c_vec",    [], ... % Psi's                                 | Vector per stage
        "R_c_vec",      [], ... % Mean radii degrees of reaction        | Vector per absolute station
        "r_hub_vec",    [], ... % Hub radii                             | Vector per absolute station
        "r_shroud_vec", [], ... % Shroud radii                          | Vector per absolute station
        "U_r_mean_vec", [], ... % U at mean radius                      | Vector per absolute station
        "r_mean_vec",   [], ... % Mean radii                            | Vector per absolute station
        "mean_index",   [], ... % Index of mean stream surface
        "num_stations", [], ... % Number of absolute stations
        "num_surfaces", [], ... % Number of stream surfaces
        "ang_vel",      [], ...`% Angular velocity
        "T0",           [], ... % Total temperature                     | Spanwise vector
        "P0",           [], ... % Total pressure                        | Spanwise vector
        "cp",           [], ... % Cp
        "gamma",        [], ... % Gamma
        "R",            [], ... % Gas constant
        "target_m_dot", [] ...  % Target mass flux
    );

    station_feeder = repmat(feeder_template, 10, 1);
end