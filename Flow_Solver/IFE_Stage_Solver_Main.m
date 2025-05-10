function stage = IFE_Stage_Solver_Main(P1, P01, T1, T01, A1, P3, gamma, R, mdot, u, react)
    % ===== RELEVANT INPUTS =====
    % P1    = Static pressure at Station 1      
    % P01   = Total pressure at Station 1
    % T1    = Static temperature at Station 1
    % T01   = Total temperature at Station 1
    % A1    = Annular area at Station 1
    % P3    = Desired pressure at exit of stage
    % Gamma = Specific heat ratio
    % R     = Gas Constant
    % mdot  = Mass flow
    % u     = Rotational velocity (scalar)
    % react = Degree of Reaction

    stage.P2 = P1*(1-(1-react)*(1-(P3/P1)^((gamma-1)/gamma))^(gamma/(gamma-1)));

    stage = IFE_Station2_Solver(stage, stage.P01, stage.T01, stage.A1, stage.P2,           stage.gamma, stage.R, stage.mdot, stage.u);
    stage = IFE_Station3_Solver(stage,                       stage.A2, stage.P2, stage.P3, stage.gamma, stage.R, stage.mdot, stage.u);
    
    
    
end

function stage = IFE_Station2_Solver(stage, P01, T01, A1, P2, gamma, R, mdot, u)
    % ===== INPUTS =====
    % stage = Struct that holds information
    % P01   = Total pressure at Station 1
    % T01   = Total temperature at Station 1
    % A1    = Annular area at Station 1
    % P2    = Desired pressure at Station 2
    % Gamma = Specific heat ratio
    % R     = Gas Constant
    % mdot  = Mass flow
    % u     = Rotational velocity

    % Assume adiabatic, exit total pressure and temp equal to input
    P02 = P01;
    T02 = T01;

    %% Determine exit mach number and exit static temperature
    % Calculate Station Y mach number (uses Equation 6)
    M2 = sqrt((2*(P02/P2)^((gamma-1)/gamma)-2)/(gamma-1));
    % Calculate Station Y static temperature (uses Equation 4)
    T2 = T02 * IFE.T_ratio_from_P_ratio(P2/P02,gamma);

    %% Determine Velocity Magnitude
    % Calculate speed of sound at Station Y
    a2 = sqrt(gamma * R * T2);
    % Calculate v_mag
    v2_mag = M2*a2;

    %% Calculate Velocity Direction
    % Density at Station Y
    rho_2 = P2/(R*T2);
    % Exit area
    A2 = mdot/(rho_2*v2_mag);
    % Alpha
    alpha2 = acosd(A2/A1);

    %% Vectors (boo-yeah)
    V = [v2_mag*cosd(alpha2), v2_mag*sind(alpha2)];
    U = [0,u];
    W = V - U;

    stage.P02 = P02;
    stage.T02 = T02;
    stage.A2  = A2;
    stage.M2  = M2;
    stage.T2  = T2;
    stage.v2_mag = v2_mag;
    stage.w2_mag = norm(W);
    stage.U2  = U;
    stage.V2  = V;
    stage.W2  = W;
end

function stage = IFE_Station3_Solver(stage, A2, P2, P3, gamma, R, mdot, u)
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

    % Converting everything into rotor rotating relative frame
    M2R = stage.w2_mag / sqrt(gamma * R * stage.T2);
    P02R = P2/( (1+(gamma-1)/2*M2R^2) ^ (-gamma/(gamma-1)) );
    T02R =  1/( (1+(gamma-1)/2*M2R^2) ^ (-gamma/(gamma-1)) );

    % Assume adiabatic, exit total pressure and temp equal to input
    P03R = P02R;
    T03R = T02R;

    T3 = T03R * IFE.T_ratio_from_P_ratio(P3/P03R, gamma);

    %% Determine exit mach number and exit static temperature
    % Calculate Station Y mach number (uses Equation 7)
    M3R = sqrt((2*(T03R/T3)-2)/(gamma-1));

    %% Determine W Velocity Magnitude
    % Calculate speed of sound at Station Y
    a3R = sqrt(gamma * R * T3);
    % Calculate w_mag
    w3_mag = M3R*a3R;

    %% Calculate Velocity Direction
    % Density at Station Y
    rho_3 = P3/(R*T3);
    % Exit area
    A3 = mdot/(rho_3*v3_mag);
    % Alpha TODO
    beta3 = acosd(A3/A2);

    %% Vectors (boo-yeah)
    W = [w3_mag*cosd(beta3), w3_mag*sind(beta3)];
    U = [0,u];
    V = W + U;

    stage.P03 = P03;
    stage.T03 = T03;
    stage.A3  = A3;
    stage.M3  = M3;
    stage.T3  = T3;
    stage.v3_mag = v3_mag;
    stage.w3_mag = norm(W);
    stage.U3  = U;
    stage.V3  = V;
    stage.W3  = W;
end