function [res, r_grid, y_grid] = quasinormal_analysis( ...
    stage_num, ...      % Stage number
    stage_123, ...      % Position within stage (local station number)
    phi_c_vec, ...      % Phi's                                 | Vector per stage
    psi_c_vec, ...      % Psi's                                 | Vector per stage
    R_c_vec, ...        % Mean radii degrees of reaction        | Vector per absolute station
    r_hub_vec, ...      % Hub radii                             | Vector per absolute station
    r_shroud_vec, ...   % Shroud radii                          | Vector per absolute station
    U_r_mean_vec, ...   % U at mean radius                      | Vector per absolute station
    r_grid, ...         % Grid of radii
    y_grid, ...         % Grid of meridional y
    r_mean_vec, ...     % Mean radii                            | Vector per absolute station
    mean_index, ...     % Index of mean stream surface
    num_stations, ...   % Number of absolute stations
    num_surfaces, ...   % Number of stream surfaces
    ang_vel, ...`       % Angular velocity
    T0, ...             % Total temperature                     | Spanwise vector
    P0, ...             % Total pressure                        | Spanwise vector
    cp, ...             % Cp
    gamma, ...          % Gamma
    R, ...              % Gas constant
    target_m_dot ...    % Target mass flux
    )

    abs_station = (stage_num-1)*2+stage_123;

    current_m_dot = target_m_dot + 1;
    new_radius_ish = r_hub_vec(abs_station);
    while abs(current_m_dot - target_m_dot) > 0.001

        if stage_123 == 1
            Wm_init = phi_c_vec(stage_num) * U_r_mean_vec(abs_station);
        elseif stage_123 == 2
            Wm_init_1 = phi_c_vec(stage_num) * U_r_mean_vec(abs_station);
            Wm_init_3 = phi_c_vec(stage_num + 1) * U_r_mean_vec(abs_station);
            Wm_init = (Wm_init_1+Wm_init_3)/2;
        elseif stage_123 == 3
            Wm_init = phi_c_vec(stage_num + 1) * U_r_mean_vec(abs_station);
        end

        r_hub_vec(abs_station:end) = new_radius_ish;
        [r_grid, y_grid, r_mean_vec, U_r_mean_vec, mean_index] = set_grid(r_shroud_vec, r_hub_vec, num_stations, num_surfaces, ang_vel);
    
        if stage_123 == 1
            C_theta = txtbk.eleven_four(R_c_vec(stage_num), psi_c_vec(stage_num), U_r_mean_vec(abs_station), r_mean_vec(abs_station), r_grid{abs_station}, 1, 1);
        elseif stage_123 == 2
            C_theta_1 = txtbk.eleven_four(R_c_vec(stage_num), psi_c_vec(stage_num), U_r_mean_vec(abs_station-1), r_mean_vec(abs_station-1), r_grid{abs_station-1}, 1, 1);
            C_theta = txtbk.eleven_five(C_theta_1, psi_c_vec(stage_num), U_r_mean_vec(abs_station-1), r_mean_vec(abs_station-1), r_grid{abs_station-1}, r_grid{abs_station});
        elseif stage_123 == 3
            C_theta = txtbk.eleven_six(R_c_vec(stage_num+1), psi_c_vec(stage_num+1), U_r_mean_vec(abs_station), r_mean_vec(abs_station), r_grid{abs_station}, 1, 1);
        end
        
        U = ang_vel .* r_grid{abs_station};
        W_theta = C_theta - U;
    
        C = sqrt(Wm_init.^2 + C_theta.^2);
        H0 = cp*T0;                              % Total enthalpy        | Assumed constant spanwise
    
        L = ones(size(U)).*0.03;            % Entropy               | Spanwise distribution        | Constant every rotation (for now)
        I = H0 - U.*C_theta;                % Rothalpy              | Spanwise distribution        | Changes every rotation  (dependent on C_theta_1
        S = zeros(size(U));                 % Entropy               | Spanwise distribution        | Constant every rotation (for now, depends on L_1??)
        T = T0 - C.^2./(2*cp);              % Static temperature    | Spanwise distribution        | Changes evert rotation  (dependent on C_!
    
        W_m = seven_fifteen(y_grid{abs_station}, r_grid{abs_station}, C_theta, W_theta, I, S, T, Wm_init, mean_index);
        [rho, rho_T, rho_P] = thermos(W_m, C_theta, T0, P0, cp, gamma, R);
        [~, new_radius_ish, current_m_dot] = annulus_adjust(y_grid{abs_station}, r_grid{abs_station}, rho, W_m, target_m_dot);
        diff = current_m_dot-target_m_dot;
        asdf = abs(diff);
        fprintf("%3f\n", asdf)
        
    end


    res = struct( ...
            "r_hub_vec", r_hub_vec, ...
            "r_shroud_vec", r_shroud_vec, ...
            "r_mean_vec", r_mean_vec, ...
            "U_r_mean_vec", U_r_mean_vec, ...
            "mean_index", mean_index, ...
            "W_m", W_m, ...
            "C_theta", C_theta, ...
            "W_theta", W_theta, ...
            "U", U, ...
            "Loss", L, ...
            "rho", rho, ...
            "rho_T", rho_T, ...
            "rho_P", rho_P, ...
            "m_dot_cur", current_m_dot ...
        );
end



function [rho, T, P] = thermos(W_m, C_theta, T0, P0, cp, gamma, R)
    C_m = W_m;
    C = sqrt(C_m.^2 + C_theta.^2);

    T = T0 - C.^2./(2*cp);                  % Static temperature    | Spanwise distribution     | Kelvin
    T_pos_test = T < 0;
    allPos = sum(T_pos_test) == 0; % If allPos is true, all of T is positive
    T_real_test = ~isreal(T); % if T_real_test has ones, then there are complex
    allReal = sum(T_real_test) == 0; % If allReal is true, then all of T are real
    if ~allPos || ~allReal
        fprintf("U too big -> ")
    end
    P = P0 .* (T./T0).^(gamma/(gamma-1));   % Static pressure       | Spanwise distribution     | Pascals
    rho = P./(R.*T);                        % Density               | Spanwise distribution     
end

function [r_new, new_radius_ish, current_m_dot] = annulus_adjust(y, r, rho, W_m, target_m_dot)
    current_m_dot = 0;
    % ==== Assumptions ====
    Kw = 1;                 % Ignoring curvature for now
    ep = zeros(size(y));    % Vertical quasi-normals

    for i = 1:length(y)-1
        dA = 2*pi*r(i)*Kw*cosd(ep(i))*(y(i+1)-y(i));
        current_m_dot = current_m_dot + rho(i) * W_m(i) * dA;
    end

    r_s = r(end);
    r_h = r(1);
    annulus_current = pi * (r_s^2 - r_h^2);
    annulus_new = annulus_current * target_m_dot / current_m_dot;

    syms r_new

    r_new_sol = solve(pi*r_s^2-pi*r_new^2 == annulus_new, r_new);
    r_new = double(r_new_sol(2));
    if isreal(r_new) == false
        if annulus_new > (pi*r_s^2)
            fprintf("shroud radius too smol -> ")
        end
        fprintf("complex radii -> ")
    end

    % new_radius = sqrt(r_s^2 - annulus_new/pi);
    new_radius_ish = r_h + (3/4) * (r_new - r_h);
end

function W_m = seven_fifteen(y, r, C_theta, W_theta, I, s, T, Wm_init, mean_index)
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

    % ==== Ong it's integration time ====
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