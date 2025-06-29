function [station_res, r_grid, y_grid, next_stop] = stage(station_info, r_grid, y_grid)
    next_stop = station_info;
    abs_station = (station_info.stage_num-1)*2+station_info.stage_123;

    %% ======== Station 1 Zoom Time ========
    [station_res, r_grid, y_grid] = quasinormal_analysis( ...
        station_info.stage_num, ...      % Stage number
        station_info.stage_123, ...      % Position within stage (local station number)
        station_info.phi_c_vec, ...      % Phi's                                 | Vector per stage
        station_info.psi_c_vec, ...      % Psi's                                 | Vector per stage
        station_info.R_c_vec, ...        % Mean radii degrees of reaction        | Vector per absolute station
        station_info.r_hub_vec, ...      % Hub radii                             | Vector per absolute station
        station_info.r_shroud_vec, ...   % Shroud radii                          | Vector per absolute station
        station_info.U_r_mean_vec, ...   % U at mean radius                      | Vector per absolute station
        r_grid, ...                      % Grid of radii
        y_grid, ...                      % Grid of meridional y
        station_info.r_mean_vec, ...     % Mean radii                            | Vector per absolute station
        station_info.mean_index, ...     % Index of mean stream surface
        station_info.num_stations, ...   % Number of absolute stations
        station_info.num_surfaces, ...   % Number of stream surfaces
        station_info.ang_vel, ...`       % Angular velocity
        station_info.T0, ...             % Total temperature                     | Spanwise vector
        station_info.P0, ...             % Total pressure                        | Spanwise vector
        station_info.cp, ...             % Cp
        station_info.gamma, ...          % Gamma
        station_info.R, ...              % Gas constant
        station_info.target_m_dot ...    % Target mass flux
        );
    
    % Plotting Station 1 Triangles
    hold on
    mean_index = station_res.mean_index;
    
    plot_triangle_span([(abs_station-1)*300,0], station_info.num_surfaces, station_res.U, station_res.W_m, station_res.W_theta, station_res.C_theta, 1)
    
    triangle_mean = {[0, station_res.U(mean_index)], [station_res.W_m(mean_index), station_res.W_theta(mean_index)], [station_res.W_m(mean_index), station_res.C_theta(mean_index)]};
    plot_vel_triangle([(abs_station-1)*300,0], triangle_mean, 1, '-g', '-b', '-r')
    
    %% ======== Linking to Next Station ========
    if station_info.stage_123 == 2
        stage_num_next = station_info.stage_num + 1;
        stage_123_next = 1;
    elseif station_info.stage_123 == 1
        stage_num_next = station_info.stage_num;
        stage_123_next = 2;
    else
        fprintf("bro i thought we agreed no station 3's")
    end
    abs_station_next = (stage_num_next-1)*2+stage_123_next;


    if station_info.stage_123 == 1
        U_curr = station_res.U;
        U_next = station_info.ang_vel .* r_grid{2};
        C_theta_next = txtbk.eleven_five(station_res.C_theta, station_info.psi_c_vec(station_info.stage_num), station_res.U_r_mean_vec(abs_station), station_res.r_mean_vec(abs_station), r_grid{abs_station}, r_grid{abs_station_next});
    elseif station_info.stage_123 == 2
        U_curr = zeros([1, station_info.num_surfaces]);
        U_next = zeros([1, station_info.num_surfaces]);
        C_theta_next = txtbk.eleven_six(station_info.R_c_vec(stage_num_next), station_info.psi_c_vec(stage_num_next), station_res.U_r_mean_vec(abs_station_next), station_res.r_mean_vec(abs_station_next), r_grid{abs_station_next}, 1, 1);
    end
    
    [T0_next, P0_next, T_next, P_next, T_cur, P_cur] = thermobridge(station_res.C_theta, station_res.W_m, U_curr, C_theta_next, station_res.W_m, U_next, station_info.T0, station_info.P0, station_info.cp, station_info.gamma, station_res.Loss);
    
    fprintf("\n")

    thermo = struct( ...
        "T0_next", T0_next, ...
        "P0_next", P0_next, ...
        "T0_cur", station_info.T0, ...
        "P0_cur", station_info.P0, ...
        "T_next",  T_next, ...
        "P_next",  P_next, ...
        "T_cur",   T_cur, ...
        "P_cur",   P_cur ...
        );

    next_stop.stage_num = stage_num_next;
    next_stop.stage_123 = stage_123_next;
    next_stop.r_hub_vec = station_res.r_hub_vec;
    next_stop.r_shroud_vec = station_res.r_shroud_vec;
    next_stop.U_r_mean_vec = station_res.U_r_mean_vec;
    next_stop.r_mean_vec = station_res.r_mean_vec;
    next_stop.mean_index = station_res.mean_index;
    next_stop.T0 = thermo.T0_next;
    next_stop.P0 = thermo.P0_next;

    station_res.thermo = thermo;
end

%% ======== Helper Functions ========
function [T02, P02, T2, P2, T1, P1] = thermobridge(C_theta_1, W_m_1, U1, C_theta_2, W_m_2, U2, T01, P01, cp, gamma, loss)

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
    % P1 = P01 - (rho .* (C1.^2)) ./ 2;
    P1 = P01 .* (T1./T01).^(gamma/(gamma-1));

    T1R = T1;
    P1R = P1;

    T01R = T1R + W1.^2./(2.*cp);
    % P01R = P1R + rho.*W1.^2./2;
    P01R = P1R .* (T1R./T01R).^(-gamma/(gamma-1)); 

    % ======== Station 2 ========
    T02Ri = T01R;                       % Assume adiabatic
    P02Ri = P01R;                       % Assume isentropic

    T02R = T02Ri;                       % Assumption becomes reality (O_o)
    P02R = P02Ri - loss.*(P01R - P1R);  % Factor in loss

    T2R = T02R - W2.^2./(2.*cp);
    % P2R = P02R - rho .* W2.^2 ./ 2;
    P2R = P02R .* (T2R./T02R).^(gamma/(gamma-1));

    T2 = T2R;
    P2 = P2R;

    T02 = T2 + C2.^2./(2.*cp);
    % P02 = P2 + rho.*C2.^2./2;
    P02 = P2 .* (T2./T02).^(-gamma/(gamma-1));

end

function plot_vel_triangle(origin, triangle, scale, U_style, W_style, C_style)
    U = triangle{1}.*scale;
    W = triangle{2}.*scale;
    C = triangle{3}.*scale;
    quiver(origin(1), origin(2), W(1), W(2), 0, W_style)
    quiver(origin(1), origin(2), C(1), C(2), 0, C_style)
    quiver(origin(1)+W(1), origin(2)+W(2), U(1), U(2), 0, U_style)
end
function plot_triangle_span(origin, num_surfaces, U_vec, W_m, W_theta, C_theta, scale)
    for i = 1:num_surfaces
        grey = (i/num_surfaces)/2 + 1/4;
        U = [0     , U_vec(i)      ].*scale;
        W = [W_m(i), W_theta(i)].*scale;
        C = [W_m(i), C_theta(i)].*scale;
        quiver(origin(1)     , origin(2)     , W(1), W(2), 0, Color=[grey, grey, grey])
        quiver(origin(1)     , origin(2)     , C(1), C(2), 0, Color=[grey, grey, grey])
        quiver(origin(1)+W(1), origin(2)+W(2), U(1), U(2), 0, Color=[grey, grey, grey])
    end
end