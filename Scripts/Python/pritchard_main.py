import math as m
import warnings
import pritchardPoints
import pritchardCurves
C:\Users\josie\Documents\PURPL\Turbojet\Turbofan\Scripts\Python\pritchard_main.py
def deg2rad(degrees):
    radians = degrees/180*m.pi
    return(radians)

def pritchard_main(params):
    failcode = "success!"
    max_iter = 300

    params["beta_IN"] = deg2rad(params["beta_IN"])
    params["beta_OUT"] = deg2rad(params["beta_OUT"])
    params["ep_IN"] = deg2rad(params["ep_IN"])
    params["ep_OUT"] = deg2rad(params["ep_OUT"])
    params["zeta"] = deg2rad(params["zeta"])

    o = 2 * m.pi * params["R"] / params["N_B"] * m.cos(params["beta_OUT"] - 2 * params["R_TE"])

    factor = 1
    ttc = params["ttc"]/100
    iter_counter = 0
    while factor > params["iteration_threshold"]:
        [pts, failcode] = pritchardPoints(params["R"], params["R_LE"], params["R_TE"], 
                                          params["Cx"], params["Ct"], params["zeta"], 
                                          params["beta_IN"], params["ep_IN"], params["beta_OUT"], 
                                          params["ep_OUT"], params["N_B"], params["o"])
        try:
            blade = pritchardCurves(pts["x_coords"], pts["y_coords"], pts["betas"], 
                                    params["R_LE"], params["R_TE"], o, params["res"])
            
        except:
            warnings.warn("bro what the hell man")
            blade["failcode"] = "knots broken"
            return
    
        blade["x_comb"], blade["y_comb"] = combiner(blade)

        if params["Ct"] ==0:
            params["Ct"] = pts["Ct"]
        
        [t_max, pos_ss, pos_ps] = max_t(blade)
        chord = m.sqrt((params["Ct"])**2 + params["Cx"]**2)
        factor = abs(t_max/chord - ttc)

        if factor > params["iteration_threshold"]:
            params["Ct"] = params["Ct"] * (3 + t_max / chord / ttc) / 4
        
        iter_counter += 1

        if iter_counter > max_iter:
            failcode = "ttc didn't converge"
            break
    
    blade["x_thicc"] = [blade["x_pressure"[pos_ps - 1]], blade["x_suction"[pos_ss -1]]]
    blade["y_thicc"] = [blade["y_pressure"[pos_ps - 1]], blade["y_suction"[pos_ss-1]]]

    #Extra params
    blade["parameters"]["R"] = params["R"]
    blade["paramters"] = params
    blade["parameters"]["o"] = o



        

            