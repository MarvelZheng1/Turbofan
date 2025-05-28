clear;clc;close all
load("Turbofan.mat")


%% ==== User-Defined Design Variables ====

rpm = 40000;

num_stages = 4;
num_stations = num_stages * 2 + 1;

num_surfaces = 50;
m_dot = 0.5;

%% ==== Inputs from the sizing script ====
T0_in  = Turbofan.Thermos.S2.T0;
T0_out  = Turbofan.Thermos.S25.T0;
P0_in = Turbofan.Thermos.S2.P0;
gamma  = Turbofan.Specs.Gammas.c_lp;
R      = Turbofan.Specs.Info.Ra;
cp     = Turbofan.Cp.c_LP;


% ==== Dimensionless Performance Parameters, per stage ====
phi_c_vec = [0.7, 0.6, 0.5, 0.4];        % Flow coefficient
psi_c_vec = [0.3, 0.3, 0.3, 0.3];        % Work coefficient ----> to be determined from temperature rises later
R_c_vec   = [0.5, 0.5, 0.5, 0.5];        % Degree of reaction


%% ==== Startup Fluff ====
% Hub and Shroud Geometry Initialization
% Vector along axis
r_shroud_vec = ones(1, num_stations)*0.001* 150;
r_hub_vec    = ones(1, num_stations)*0.001* 75;

% Mean Radius Initialization
r_mean_vec   = ones(1, num_stations);
r_mean_real_vec   = ones(1, num_stations);

% Misc
ang_vel = rpm * 2 * pi / 60;
% y = r - r_hub
r_grid = cell(1,num_stations); % {1} = leftmost quasi-normal  |  (1) = hub radius  |  Basically {x}(y) is cartesian
y_grid = cell(1,num_stations); % {1} = leftmost quasi-normal  |  (1) = hub radius  |  Basically {x}(y) is cartesian

% Find meanline radii
syms r_mean_sym;
for i = 1:num_stations
    r_mean_sols = solve(r_shroud_vec(i)^2-r_mean_sym^2 == r_mean_sym^2-r_hub_vec(i)^2, r_mean_sym);
    r_mean_real_vec(i) = double(r_mean_sols(2));
    r_grid{i} = linspace(r_hub_vec(i), r_shroud_vec(i), num_surfaces);         % Grid of y-values, to be used in flow-field analysis
end

for i = 1:num_stations
    y_grid{i} = r_grid{i} - r_grid{i}(1);
end

r_mean_dist = r_grid{1}-r_mean_real_vec(1);
[~, mean_index] = min(abs(r_mean_dist));
for i = 1:num_stations
    r_mean_vec(i) = r_grid{i}(mean_index);                       % Y-value streamline closest to the mean radius of station 1 (equal annulus area above and below mean radius)
end


clear i r_mean_dist r_mean_sols r_mean_sym

% Find meanline U velocities (returns mean radius values of U axially)
U_r_mean_vec = ang_vel.*r_mean_vec;


% ==== Wm_init Initialization ====
Wm_init_1 = phi_c_vec(1) * U_r_mean_vec(1);
Wm_init_2 = phi_c_vec(2) * U_r_mean_vec(3);
Wm_init_3 = (Wm_init_1 + Wm_init_2) / 2;

C_theta_1 = eleven_four(R_c_vec(1), psi_c_vec(1), U_r_mean_vec(1), r_mean_vec(1), r_grid{1}, 1, 1);
U_1 = ang_vel .* r_grid{1};
W_theta_1 = C_theta_1 - U_1;

% ==== First-Pass LIST Initialization ====
C_init = sqrt(Wm_init_1.^2 + C_theta_1.^2);
h0 = cp*T0_in;                              % Total enthalpy        | Assumed constant spanwise

L_init = ones(size(U_1)).*0.03;             % Entropy               | Spanwise distribution
I_init = h0 - U_1.*C_theta_1;               % Rothalpy              | Spanwise distribution
S_init = zeros(size(U_1));                  % Entropy               | Spanwise distribution
T_init = T0_in - C_init.^2./(2*cp);         % Static temperature    | Spanwise distribution


W_m_1 = seven_fifteen(y_grid{1}, I_init, S_init, r_grid{1}, T_init, C_theta_1, W_theta_1, Wm_init_1, mean_index);
stat1_thermos = station_thermodynamics(W_m_1, C_theta_1, U_1, S_init, T0_in, P0_in, cp, gamma, L_init, R);
adjusted_hub_radii = annulus_adjust(y_grid{1}, r_grid{1}, stat1_thermos.rho, W_m_1);



function thermos = station_thermodynamics(W_m, C_theta, U, S, T0, P0, cp, gamma, loss, R)
    C_m = W_m;
    C = sqrt(C_m.^2 + C_theta.^2);

    h0 = cp*T0;                             % Total enthalpy        | Assumed constant spanwise
    T = T0 - C.^2./(2*cp);                  % Static temperature    | Spanwise distribution
    P = P0 .* (T./T0).^(gamma/(gamma-1));   % Static pressure       | Spanwise distribution
    rho = P./(R.*T);                        % Density               | Spanwise distribution
    dS = loss.*0.5.*U.^2./T0;               % Entropy change        | Spanwise distribution
    S = S + dS;                             % Entropy               | Spanwise distribution
    I = h0 - U.*C_theta;                     % Rothalpy              | Spanwise distribution

    % ==== Delicious Thermos Meal Recipie ====
    thermos.h0 = h0;
    thermos.T = T;
    thermos.P = P;
    thermos.rho = rho;
    thermos.dS = dS;
    thermos.S = S;
    thermos.I = I;
end

function current_m_dot = annulus_adjust(y, r, rho, W_m)
    current_m_dot = 0;
    % ==== Assumptions ====
    Kw = 1;                 % Ignoring curvature for now
    ep = zeros(size(y));    % Vertical quasi-normals

    for i = 1:length(y)-1
        dA = 2*pi*r(i)*Kw*cosd(ep(i))*(y(i+1)-y(i));
        current_m_dot = current_m_dot + rho(i) * W_m(i) * dA;
    end
end

function C_theta_1 = eleven_four(Rc_c1, psi_c1, U_c1, r_c1, r_1, n, m)
    C_theta_1 = U_c1.*((1-Rc_c1).*(r_c1./r_1).^n - psi_c1./2.*(r_c1./r_1).^m);
end

function W_m = seven_fifteen(y, I, s, r, T, C_theta, W_theta, Wm_init, mean_index)
    % ==== Input Variables ====
    %       y | spanwise coordinate (vector)
    %       I | rothalpy distribution (vector)
    %       s | entropy distribution (vector)
    %       r | radii (vector)
    %       T | total temperature (scalar)
    % C_theta | absolute tangential velocity (vector)
    % W_theta | relative tangential velocity (vector)
    % Wm_init | known boundary value of Wm at starting y (hub)

    % ==== Assumptions ====
    Km = 0;                     % No streamline curvature
    dWm_dm = zeros(size(y));    % No streamline variation
    ep = zeros(size(y));        % Vertical quasi-normals

    % ==== Gradients ====
    dI_dy = gradient(I, y);
    ds_dy = gradient(s, y);
    drC_dy = gradient(r.*C_theta, y);

    % ==== Used later when we want to consider streamline variation ====
    % dphi_dy = gradient(phi, y);
    % dWm_dm = Wm./(1-Mm.^2) .* (-(1+M_theta.^2) .* sind(phi)./r - 1/cosd(ep).*dphi_dy - Km.*tand(ep));

    % ==== Initialization ====
    W_m = ones(size(y))*Wm_init;

    % ==== Le f's ====
    f1 = -Km.*cosd(ep) + sind(ep)./W_m .* dWm_dm;   % zero for initial sizing
    f3 = dI_dy - T.*ds_dy - (W_theta./r).*drC_dy;

    % ==== Ong it's marching time
    for i = mean_index:length(y)-1
        dy = y(i+1)-y(i);
        dWm_dy = f1(i)*W_m(i) + f3(i)/W_m(i);
        W_m(i+1) = W_m(i) + dy * dWm_dy;
    end

    for i = mean_index:-1:2
        dy = y(i)-y(i-1);
        dWm_dy = f1(i)*W_m(i) + f3(i)/W_m(i);
        W_m(i-1) = W_m(i) - dy * dWm_dy;
    end
end