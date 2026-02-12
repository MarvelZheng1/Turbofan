import math as m
import warnings
import pritchardPoints
import pritchardCurves
from REF_structs import params, parameters, blade, pts

def deg2rad(degrees):
    radians = degrees/180*m.pi
    return(radians)

def rad2deg(radians):
    degrees = radians*180/m.pi
    return(degrees)

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
    blade["parameters"] = params
    blade["parameters"]["o"] = o
    blade["parameters"]["pitch"] = 2*m.pi*params["R"]/params["N_B"]
    blade["parameters"]["t_max"] = max_t(blade)
    blade["parameters"]["t_min"] = min_t(blade)
    if blade["parameters"]["t_min"] < 1:
        failcode = "blade too thin"
    
    blade["parameters"]["zweifel"] = (4*m.pi*params["R"]) / (params["Cx"]*params["N_B"]) * m.sin(params["beta_IN"] - params["beta_OUT"]) * m.cos(params["beta_OUT"]) / m.cos(params["beta_IN"])
    blade["parameters"]["blockage_IN"] = 2 * params["R_LE"] / blade["params"]["pitch"] * m.cos(params["beta_IN"]) * 100
    blade["parameters"]["blockage_OUT"] = 2*params["R_TE"] / (blade["parameters"]["pitch"] * m.cos(params["beta_OUT"])) * 100 
    blade["parameters"]["chord"] = m.sqrt(params["Ct"]^2 + params["Cx"]^2)
    blade["parameters"]["calc_ttc"] = blade["parameters"]["t_max"] / blade["parameters"]["chord"]

    """""
    Notes from 2/11 for 2/12:
    -Continued on conversion to python with dictionaries
    -Wanted to wait to figure out dataclasses until in person


    """
    

    #unconverted matlab copy paste
    blade.parameters.pitch_to_chord = blade.parameters.pitch/blade.parameters.chord;    
    blade.parameters.height_to_chord = blade.parameters.blade_height/blade.parameters.chord;

    blade.parameters.beta_IN = rad2deg(blade.parameters.beta_IN);
    blade.parameters.beta_OUT = rad2deg(blade.parameters.beta_OUT);
    blade.parameters.ep_IN = rad2deg(blade.parameters.ep_IN);
    blade.parameters.ep_OUT = rad2deg(blade.parameters.ep_OUT);
    blade.parameters.zeta = rad2deg(blade.parameters.zeta); 
