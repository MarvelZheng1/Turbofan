from dataclasses import dataclass
import numpy as np
import math as m


import REF_AEQ
import REF_structs
import Station_Thermo
import Component_Sizing


# Dataclasses for storing important component parameters and values

# Initializing dataclasses for component efficiencies, specific heat ratios, and pressure ratios
eta = REF_structs.ByComponent(None,
                  0.94, # Diffuser
                  0.85, # Fan
                  0.98, # Fan Nozzle
                  0.83, # LP Compressor
                  0.83, # HP Compressor
                  1.00, # Burner
                  0.89, # HP Turbime
                  0.89, # LP Turbine
                  0.98  # Nozzle
)
gamma = REF_structs.ByComponent(1.4,  # Ambient
                    1.4,  # Diffuser
                    1.4,  # Fan
                    1.4,  # Fan Nozzle
                    1.4,  # LP Compressor
                    1.4,  # HP Compressor
                    1.3,  # Burner
                    1.32, # HP Turbine
                    1.32, # LP Turbine
                    1.34  # Nozzle
)
Pr = REF_structs.ByComponent(None,
                 None,
                 1.2,        # Fan
                 None,
                 m.sqrt(5), # LP Compressor
                 m.sqrt(5), # HP Compressor
                 1,          # Burner
                 None,
                 None,
                 None
)

# ======== Ambient Conditions and General Engine Parameters ========
cycle_params = {
    "eta"    : eta,             # nondimensional
    "gamma"  : gamma,           # nondimensional
    "Pr"     : Pr,              # nondimensional
    "T_0"    : 298,             # TODO | At ambient, v=0, so T0 = T
    "P_0"    : 101300,          # TODO | At ambient, v=0, so P0 = P
    "M_f"    : 0.85,            # nondimensional
    "Ra"     : 287,             # TODO
    "Rp"     : 287,             # TODO
    "QR"     : 45000000,        # TODO
    "bypass" : 3,               # nondimensional
    "combustion_temp": 1750,    # TODO
}
thermo, Cps, thrust, eta_calc = Station_Thermo.thermoCalcs(cycle_params)

compressor_params = {
    "gamma"       : gamma.cLP,
    "Cp_cLP"      : Cps.cLP,
    "T0_1"        : thermo.S2.T0,
    "P0_1"        : thermo.S2.P0,
    "Pr"          : Pr.cLP,
    "e_c"         : 0.9,            # Polytropic Efficiency
    "httrr"       : 0.6,            # Hub-to-tip radius ratio
    "deHaller"    : 0.72,           # De Haller's Criterion value
    "min_Re"      : 300000,         # For blade sizing
    "mu_kin"      : 1.46e-5,        # Kinematic viscosity
    "alpha_1m"    : 30,             # Inlet absolute flow angle, in degrees
    "Mz_1m"       : 0.45,           # Inlet axial mach number
    "U_tip_inlet" : 350,            # Inlet rotor tip speed, m/s
    "r_tip_inlet" : 100/1000        # Inlet rotor tip radius, in mm (divided by 1000 to get meters)
}
# compressor_info = Component_Sizing.Compressor_Sizing(compressor_params)
# turbine_info = Component_Sizing.Turbine_Sizing(turbine_params)

# plottingFuncs.plot_T0P0_vs_Stations(thermo)

with open("results.txt", "w") as txt:
    txt.write("======== Thrust and Efficiencies ========\n")
    txt.write("    ST |{:8.3f} N\n".format(thrust[0]))
    txt.write("  TSFC |{:8.3f} N\n".format(thrust[1]))
    txt.write(" eta_p |{:8.3f} %\n".format(eta_calc[0]*100))
    txt.write("eta_th |{:8.3f} %\n".format(eta_calc[1]*100))
    txt.write(" eta_0 |{:8.3f} %\n".format(eta_calc[2]*100))