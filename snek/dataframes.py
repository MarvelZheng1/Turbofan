from dataclasses import dataclass

@dataclass
class fullVelTriInfo:
    C_1m: float
    C_2m: float
    C_3m: float
    W_1m: float
    W_2m: float
    W_3m: float
    U_1m: float
    U_2m: float
    U_3m: float
    z_1m: float
    z_2m: float
    z_3m: float
    
    Mc_1m: float
    Mc_2m: float
    Mc_3m: float
    Mw_1m: float
    Mw_2m: float
    Mw_3m: float
    Mz_1m: float
    Mz_2m: float
    Mz_3m: float
    
    Ctheta_1m: float
    Ctheta_2m: float
    Ctheta_3m: float
    Wtheta_1m: float
    Wtheta_2m: float
    Wtheta_3m: float
    
    alpha_1m: float
    alpha_2m: float
    alpha_3m: float
    beta_1m: float
    beta_2m: float
    beta_3m: float