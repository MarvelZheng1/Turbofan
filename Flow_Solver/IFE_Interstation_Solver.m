function stage = IFE_Interstation_Solver(stage, type, P0x, T0x, Ax, Py, gamma, R, mdot, u)
    % ===== INPUTS =====
    % stage = Struct that holds information
    % type  = "Stator" or "Rotor"
    % Px    = Static pressure at Station X      
    % P0x   = Total pressure at Station X
    % Tx    = Static temperature at Station X
    % T0x   = Total temperature at Station X
    % Ax    = Annular area at Station X
    % Py    = Desired pressure at exit
    % Gamma = Specific heat ratio
    % R     = Gas Constant
    % mdot  = Mass flow
    % u     = Rotational velocity (scalar)

    % =====OUTPUTS =====
    % v_mag = Absolute velocity magnitude
    % w_mag = Relative velocity magnitude
    % V     = Absolute velocity vector
    % W     = Relative velocity vector

    % Assume adiabatic, exit total pressure and temp equal to input
    P0y = P0x;
    T0y = T0x;

    %% Determine exit mach number and exit static temperature
    % Calculate Station Y mach number (uses Equation 6)
    if type == "stator"
        My = sqrt((2*(P0y/Py)^((gamma-1)/gamma)-2)/(gamma-1));
    elseif type == "rotor"
        My = stage.w2_mag / sqrt(gamma * R * stage.T2);
    % Calculate Station Y static temperature (uses Equation 4)
    Ty = T0y * (Py/P0y)^((gamma-1)/gamma);

    %% Determine Velocity Magnitude
    % Calculate speed of sound at Station Y
    a = sqrt(gamma * R * Ty);
    % Calculate v_mag
    v_mag = My*a;

    %% Calculate Velocity Direction
    % Density at Station Y
    rho_y = Py/(R*Ty);
    % Exit area
    Ay = mdot/(rho_y*v_mag);
    % Alpha
    alpha = acosd(Ay/Ax);

    %% Vectors (boo-yeah)
    V = [v_mag*cosd(alpha), v_mag*sind(alpha)];
    U = [0,u];
    W = V - U;

    if type == "stator"
        stage.P02 = P0y;
        stage.T02 = T0y;
        stage.A2  = Ay;
        stage.M2  = My;
        stage.T2  = Ty;
        stage.v2_mag = v_mag;
        stage.w2_mag = norm(W);
        stage.V2  = V;
        stage.W2  = W;
    elseif type == "rotor"
        stage.P03 = P0y;
        stage.T03 = T0y;
        stage.A3  = Ay;
        stage.M3  = My;
        stage.T3  = Ty;
        stage.V3  = V;
        stage.W3  = W;
    end


end