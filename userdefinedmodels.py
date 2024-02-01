from pvlib import modelchain
import numpy as np
import pandas as pd

def mppt(model_chain):
    """
    Calculate the MPPT output power, as a Series [W] and store in 
    model_chain.results.ac. Uses `v_batt`, `v_start`, and `v_continue` parameters 
    from `inverter_parameters` of the model_chain.system (the PVSystem)
    """
    #helper function
    def calculate_power(v_mp, v_oc, p_mp):
        power_ac = np.empty_like(p_mp)
        v_batt = model_chain.system.inverter_parameters['v_batt']
        v_start = v_batt + model_chain.system.inverter_parameters['v_start']
        v_continue = v_batt + model_chain.system.inverter_parameters['v_continue']

        on = False
        for i, (vmp, voc, pmp) in enumerate(zip(v_mp, v_oc, p_mp)):
            if on:
                if vmp < v_continue:
                    on = False
                    power_ac[i] = 0
                else:
                    power_ac[i] = pmp
            else:
                if voc >= v_start:
                    on = True
                    power_ac[i] = pmp
                
        if isinstance(p_mp, pd.Series):
            power_ac = pd.Series(power_ac, index=p_mp.index)
        return power_ac


    if isinstance(model_chain.results.dc, tuple):
        model_chain.results.ac = tuple(calculate_power(
            df['v_mp'], 
            df['v_oc'], 
            df['p_mp']
            ) for df in model_chain.results.dc)
    else:
        model_chain.results.ac = calculate_power(model_chain.results.dc['v_mp'], 
                                                 model_chain.results.dc['v_oc'], 
                                                 model_chain.results.dc['p_mp'])

    return model_chain

def pv_wire_loss(model_chain):
    return model_chain
    