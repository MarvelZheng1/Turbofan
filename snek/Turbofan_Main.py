from dataclasses import dataclass
from Station_Thermo import thermoCalcs
import numpy as np
import plottingFuncs

# Dataclasses for storing important component parameters and values
@dataclass
class ByComponent:
    a:   float      # Ambient
    d:   float      # Diffuser
    f:   float      # Fan
    fn:  float      # Fan Nozzle
    cLP: float      # LP Compressor
    cHP: float      # HP Compressor
    b:   float      # Burner
    tHP: float      # HP Turbine
    tLP: float      # LP Turbine
    n:   float      # Nozzle

# Initializing dataclasses for component efficiencies, specific heat ratios, and pressure ratios
eta = ByComponent(None,
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
gamma = ByComponent(1.4,  # Ambient
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
Pr = ByComponent(None,
                 None,
                 1.2,        # Fan
                 None,
                 np.sqrt(5), # LP Compressor
                 np.sqrt(5), # HP Compressor
                 1,          # Burner
                 None,
                 None,
                 None
)

# ======== Ambient Conditions and General Engine Parameters ========
params = {
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

thermo, thrust, eta_calc = thermoCalcs(params)
print("\n====== Thrust and Efficiencies ======")
print("    ST |{:8.3f} N".format(thrust[0]))
print("  TSFC |{:8.3f} N".format(thrust[1]))
print(" eta_p |{:8.3f} %".format(eta_calc[0]*100))
print("eta_th |{:8.3f} %".format(eta_calc[1]*100))
print(" eta_0 |{:8.3f} %".format(eta_calc[2]*100))

# plottingFuncs.plot_T0P0_vs_Stations(thermo)