import numpy as np
from dataclasses import dataclass

@dataclass
class StationTnP:
    T0: float
    P0: float
@dataclass
class StationThermo:
    Sa:  StationTnP
    S15: StationTnP
    S2:  StationTnP
    S25: StationTnP
    S3:  StationTnP
    S4:  StationTnP
    S45: StationTnP
    S5:  StationTnP
    S6:  StationTnP
    S7:  StationTnP
    S8:  StationTnP

def thermoCalcs(params):
    eta    = params["eta"]                      # Efficiencies
    gamma  = params["gamma"]                    # Specific heat ratios
    Pr     = params["Pr"]                       # Design pressure ratios
    T_0    = params["T_0"]                      # Ambient temp
    P_0    = params["P_0"]                      # Ambient pressure
    M_f    = params["M_f"]                      # Flight mach number
    Ra     = params["Ra"]                       # Gas constant of air
    Rp     = params["Rp"]                       # Gas constant of combustion products
    QR     = params["QR"]                       # Heat of reaction for combustion
    bypass = params["bypass"]                   # Bypass ratio
    combustion_temp = params["combustion_temp"] # Combustion temperature

    # lmao who cares about station 1 am i right (???)
    # ======== Station 1.5: Diffuser Outlet /Fan Inlet ========
    T0_15 = T_0*(1 + (gamma.a-1)/2 * M_f**2)
    P0_15 = P_0*(1 + eta.d*(T0_15/T_0 - 1))**(gamma.d/(gamma.d-1))

    # ======== Station 2: Fan Outlet/LP Compressor Inlet ========
    Cp_f = (Ra*gamma.f)/(gamma.f-1)     # Specific heat of fan

    T0_2 = T0_15*(1 + 1/eta.f*(Pr.f**((gamma.f-1)/gamma.f)-1))
    P0_2 = P0_15*Pr.f

    # ======== Station 2.5: LP Compressor Outlet/HP Compressor Inlet ========
    Cp_cLP = (Ra*gamma.cLP)/(gamma.cLP-1)     # Specific heat of compressor

    T0_25 = T0_2*(1 + 1/eta.cLP*(Pr.cLP**((gamma.cLP-1)/gamma.cLP)-1))
    P0_25 = P0_2*Pr.cLP

    # ======== Station 3: HP Compressor Outlet/Burner Inlet ========
    Cp_cHP = (Ra*gamma.cHP)/(gamma.cHP-1)     # Specific heat of compressor

    T0_3 = T0_25*(1 + 1/eta.cHP*(Pr.cHP**((gamma.cHP-1)/gamma.cHP)-1))
    P0_3 = P0_25*Pr.cHP

    # ======== Station 4: Burner Outlet/HP Turbine Inlet ========
    Cp_b = (Rp*gamma.b)/(gamma.b-1)     # Specific heat of burner

    T0_4 = combustion_temp
    P0_4 = P0_3*Pr.b

    fr = (T0_4/T0_3 - 1)/((eta.b*QR)/(Cp_b*T0_3)-T0_4/T0_3)   # Fuel-air ratio

    # ======== Station 4.5: HP Turbine Outlet/LP Turbine Inlet ========
    Cp_tHP = (Rp*gamma.tHP)/(gamma.tHP-1)     # Specific heat of turbine

    T0_45 = ((1+fr)*T0_4*Cp_tHP - Cp_cHP*(T0_3-T0_25)) / ((1+fr)*Cp_tHP)
    P0_45 = P0_4*(1 - 1/eta.tHP*(1 - T0_45/T0_4))**(gamma.tHP/(gamma.tHP-1))

    # ======== Station 5: LP Turbine Outlet/Nozzle Inlet ========
    Cp_tLP = (Rp*gamma.tLP)/(gamma.tLP-1)     # Specific heat of turbine

    T0_5 = ((1+fr)*T0_45*Cp_tLP - Cp_cLP*(T0_25-T0_2) - bypass*Cp_f*(T0_2-T0_15)) / ((1+fr)*Cp_tLP)
    P0_5 = P0_45*(1 - 1/eta.tLP*(1 - T0_5/T0_45))**(gamma.tLP/(gamma.tLP-1))

    # ======== Station 6: Afterburner (there is none lmao) ========
    T0_6 = T0_5
    P0_6 = P0_5

    # ======== Station 7: Nozzle Outlet ========
    T0_7 = T0_6
    P0_7 = P0_6

    # ======== Station 8: Aft Ambient ========
    T_8 = T_0
    P_8 = P_0

    # ======== Nozzle Exit Velocities ========
    u_ec = np.sqrt(2*eta.n *(gamma.n /(gamma.n -1))*Rp*T0_7*(1 - (P_8/P0_7)**((gamma.n -1)/gamma.n )))
    u_ef = np.sqrt(2*eta.fn*(gamma.fn/(gamma.fn-1))*Ra*T0_2*(1 - (P_8/P0_2)**((gamma.fn-1)/gamma.fn)))

    # ======== Performance Metrics ========
    u_a = M_f * np.sqrt(gamma.a*Ra*T_0)
    
    ST = (1+fr)*u_ec + bypass*u_ef - (1+bypass)*u_a     # Specific Thrust
    TSFC = fr/ST                                        # Thrust Specific Fuel Consumption
    eta_p  = ST*u_a / ((1+fr)*((u_ec**2)/2) + bypass*((u_ef**2)/2) - (1+bypass)*((u_a**2)/2))     # Propulsive Efficiency
    eta_th = ((1+fr)*((u_ec**2)/2) + bypass*((u_ef**2)/2) - (1+bypass)*((u_a**2)/2)) / (fr*QR)    # Thermal Efficiency
    eta_0  = eta_p*eta_th    


    T0P0 = StationThermo(StationTnP(T_0,   P_0),
                        StationTnP(T0_15, P0_15),
                        StationTnP(T0_2,  P0_2),
                        StationTnP(T0_25, P0_25),
                        StationTnP(T0_3,  P0_3),
                        StationTnP(T0_4,  P0_4),
                        StationTnP(T0_45, P0_45),
                        StationTnP(T0_5,  P0_5),
                        StationTnP(T0_6,  P0_6),
                        StationTnP(T0_7,  P0_7),
                        StationTnP(T_8,  P_8),
                        )
    
    return[T0P0, [ST, TSFC], [eta_p, eta_th, eta_0]]