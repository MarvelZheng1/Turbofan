clear;clc;close all

syms m_sym          % Meridional coordinate
syms y_sym          % Quasi-normal coordinate
syms z_sym
syms r_sym

syms Wm_sym         % Meridional velocity
syms I_sym          % Rothalpy
syms s_sym          % Entropy

%% User-Defined Design Variables:

rpm = 40000;

num_stages = 4;
num_stations = num_stages * 2 + 1;

num_surfaces = 50;

% Dimensionless Performance Parameters
phi_c_vec = [0.7, 0.7, 0.6, 0.6, 0.5, 0.5, 0.4, 0.4, 0.3];        % Flow coefficient
psi_c_vec = ones(1, num_stations)*0.3;        % Work coefficient
R_c_vec   = ones(1, num_stations)*0.5;        % Degree of reaction


%% Startup Fluff
% Hub and Shroud Geometry Initialization
% Vector along axis
r_shroud_vec = ones(1, num_stations)*150;
r_hub_vec    = ones(1, num_stations)*75;

% Mean Radius Initialization
r_mean_vec   = ones(1, num_stations);
r_mean_real_vec   = ones(1, num_stations);

% Misc
ang_vel = rpm * 2 * pi / 60;
y_grid = cell(1,num_stations);

% Find meanline radii
syms r_mean_sym;
for i = 1:num_stations
    r_mean_sols = solve(r_shroud_vec(i)^2-r_mean_sym^2 == r_mean_sym^2-r_hub_vec(i)^2, r_mean_sym);
    r_mean_real_vec(i) = double(r_mean_sols(2));
    y_grid{i} = linspace(r_hub_vec(i), r_shroud_vec(i), num_surfaces);         % Grid of y-values, to be used in flow-field analysis
end

r_mean_dist = y_grid{1}-r_mean_real_vec(1);
[~, closest] = min(abs(r_mean_dist));
for i = 1:num_stations
    r_mean_vec(i) = y_grid{i}(closest);                       % Y-value streamline closest to the mean radius of station 1 (equal annulus area above and below mean radius)
end

clear i r_mean_dist r_mean_sols r_mean_sym

% Find meanline U velocities (returns mean radius values of U axially)
U_r_mean_vec = ang_vel.*r_mean_vec;

% Find 
[Stat1, Stat2, Stat3] = velocity_triangle_generation(1, r_mean_vec, U_r_mean_vec, closest, phi_c_vec, psi_c_vec, R_c_vec, y_grid, false, false);






function [Loc_Stat1, Loc_Stat2, Loc_Stat3] = velocity_triangle_generation(stage_num, r_mean_vec, U_r_mean_vec, closest, phi_c_vec, psi_c_vec, R_c_vec, y_grid, plot, last)
    
    % Subscripts 1, 2, and 3 refer to stage-local station numbers:
    % 1 = rotor inlet, 2 = rotor outlet/stator inlet, 3 = stator outlet

    % Setting local, stage station numbers
    station_r_in = 2*(stage_num-1) + 1;
    station_r_s  = station_r_in + 1;
    station_s_out  = station_r_s + 1;
    
    % Reference mean radii
    r_c1 = r_mean_vec(station_r_in);
    r_c3 = r_mean_vec(station_s_out);

    % Reference mean reaction, psi (work coeff)
    Rc_c1 = R_c_vec(station_r_in);
    Rc_c3 = R_c_vec(station_s_out);
    psi_c1 = psi_c_vec(station_r_in);
    psi_c3 = psi_c_vec(station_s_out);

    % Reference U
    U_c1 = U_r_mean_vec(station_r_in);
    U_c2 = U_r_mean_vec(station_r_s);
    U_c3 = U_r_mean_vec(station_s_out);

    % Reference meridional C
    C_mc1 = U_c1 * phi_c_vec(station_r_in);
    C_mc3 = U_c3 * phi_c_vec(station_s_out);
    C_mc2 = (C_mc1 + C_mc3)/2;
    
    % Discrete radii
    r_1 = y_grid{station_r_in};
    r_2 = y_grid{station_r_s};
    r_3 = y_grid{station_s_out};
    
    % n and m for radial equilibrium (1 and 1 is free vortex)
    n = 1;
    m = 1;
    
    % Radial equilibrium'd tangential component of C
    C_theta_1 = U_c1.*((1-Rc_c1).*(r_c1./r_1).^n - psi_c1./2.*(r_c1./r_1).^m);
    C_theta_2 = (r_1.*C_theta_1 + psi_c1.*U_c1.*r_c1) ./ r_2;
    C_theta_3 = U_c3.*((1-Rc_c3).*(r_c3./r_3).^n - psi_c3./2.*(r_c3./r_3).^m);

    % constant = 1;
    % if last
    %     C_theta_3 = C_theta_2.*constant;
    % elseif
    %     C_theta_3 = U_c3.*((1-Rc_c3).*(r_c3./r_3).^n - psi_c3./2.*(r_c3./r_3).^m);
    % end

    % Velocity Triangles for stage stations 1 2 and 3
    % U1 = [0,U_c1];
    % C1 = [C_mc1, C_theta_1(closest)];
    % W1 = C1-U1;
    % 
    % U2 = [0,U_c2];
    % C2 = [C_mc2, C_theta_2(closest)];
    % W2 = C2-U2;
    % 
    % U3 = [0,U_c3];
    % C3 = [C_mc3, C_theta_3(closest)];
    % W3 = C3-U3;
    % 
    % % Function return packaging
    % Loc_Stat1 = {U1, W1, C1};
    % Loc_Stat2 = {U2, W2, C2};
    % Loc_Stat3 = {U3, W3, C3};
    % 
    % if plot
    %     figure(Name='Stage Triangle Series')
    %     hold on
    %     scale = 1e-5;
    %     plot_vel_triangle([0,0], Loc_Stat1, scale)
    %     plot_vel_triangle([4,0], Loc_Stat2, scale)
    %     plot_vel_triangle([8,0], Loc_Stat3, scale)
    % 
    %     figure(Name='Stage Triangle Superimposed')
    %     hold on
    %     scale = 1e-5;
    %     plot_vel_triangle([0,0], Loc_Stat1, scale)
    %     plot_vel_triangle([0,0], Loc_Stat2, scale)
    %     plot_vel_triangle([0,0], Loc_Stat3, scale)
    % end

    % for i = 1:num_surfaces
    %     U = [0,U_c1];
    %     C = [C_mc1, C_theta_1(i)];
    %     W = C-U;
    % 
    %     quiver(0,0, W(1), W(2), 0, '-r');
    %     quiver(0,0, C(1), C(2), 0, '-b');
    %     quiver(W(1),W(2), U(1), U(2), 0, '-k');
    % end
end

function plot_vel_triangle(origin, triangle, scale)
    U = triangle{1}.*scale;
    W = triangle{2}.*scale;
    C = triangle{3}.*scale;
    quiver(origin(1), origin(2), W(1), W(2), 0, '-b')
    quiver(origin(1), origin(2), C(1), C(2), 0, '-r')
    quiver(origin(1)+W(1), origin(2)+W(2), U(1), U(2), 0, '-k')
end


function W_m = seven_fifteen_i_guess(y, I, s, r, C_theta, phi, Wm, Mm, M_theta, ep, Km, W_theta)
    dI_dy = gradient(I, y);
    ds_dy = gradient(s, y);
    drC_dy = gradient(r.*C_theta, y);

    dphi_dy = gradient(phi, y);
    dWm_dm = Wm./(1-Mm.^2) .* (-(1+M_theta.^2) .* sind(phi)./r - 1/cosd(ep).*dphi_dy - Km.*tand(ep));

    f1 = -Km.*cosd(ep) + sind(ep)./Wm * dWm_dm;
    f3 = dI_dy - T.*ds_dy - (W_theta./r).*drC_dy;

    dA = 2.*pi.*r.*Kw.*cosd(ep).*dy;
    dm_dot = rho.*Wm.*dA;
    f4 = f3.*p.*dA./dm_dot;

    F = exp(trapz(f1));

    W_m = Constant.*F + F.*trapz(f4./F);
end








% [  ~  , Stat4, Stat5] = velocity_triangle_generation(2, r_mean_vec, U_r_mean_vec, closest, phi_c_vec, psi_c_vec, R_c_vec, y_grid, false, false);
% [  ~  , Stat6, Stat7] = velocity_triangle_generation(3, r_mean_vec, U_r_mean_vec, closest, phi_c_vec, psi_c_vec, R_c_vec, y_grid, false, false);
% [  ~  , Stat8, Stat9] = velocity_triangle_generation(4, r_mean_vec, U_r_mean_vec, closest, phi_c_vec, psi_c_vec, R_c_vec, y_grid, false, true);
% 
% Station_Vel_Triangles = {Stat1, Stat2, Stat3, Stat4, Stat5, Stat6, Stat7, Stat8, Stat9};


% scale = 1e-5;
% hold on
% for i = 1:length(Station_Vel_Triangles)
%     plot_vel_triangle([0,0], Station_Vel_Triangles{i}, scale)
% end

% figure()
% scale = 1e-5;
% hold on
% for i = 1:length(Station_Vel_Triangles)
% 
%     plot_vel_triangle([4*(i-1),0], Station_Vel_Triangles{i}, scale)
% end





















figure(Name="Streamlines")
hold on
for i = 1:num_surfaces
    streamline = ones(1,num_stations);
    for j = 1: num_stations
        streamline(j) = y_grid{j}(i);
    end
    plot(streamline, 'k--')
end
plot([r_shroud_vec',r_hub_vec', r_mean_vec'], '.-')
ylim([min(r_hub_vec)-10, max(r_shroud_vec)+10])
xlim([0, num_stations+1])

function [f1_sol, f2_sol, f3_sol] = f123(y)
    f1_sol = 0; %-Km.*cosd(ep) + (sind(ep)./Wm).*(diff(Wm_sym, m_sym));
    f2_sol = 0;
    f3_sol = 0;
end