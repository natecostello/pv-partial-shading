from pvlib import modelchain
import numpy as np
import pandas as pd

def mppt(model_chain):
    """
    Calculates the Maximum Power Point Tracking (MPPT) output power and stores it
    in `model_chain.results.ac` as a pandas Series measured in watts (W).

    Parameters:
        model_chain (ModelChain): An instance of a ModelChain that includes the
        PV system configuration and solar position. It must contain the DC model output
        and inverter parameters.

    Optional Parameters:
        - `mppt_eff` (float): Efficiency of MPPT, default 1.0 if not provided.
        - `r_batt_wire` (float): Resistance of battery wire, default 0 if not provided.
        - `v_batt` (float): Battery voltage, default 14.4 if not provided. 
        - `v_start_delta` (float): Start voltage for MPPT, default 5.0 if not provided.
        - `v_continue_delta` (float): Continue voltage for MPPT, default 1.0 if not provided.

    The function expects the use of DC models (e.g., "cec" or "desoto") that provide
    `v_mp`, `v_oc`, and `p_mp` values in `model_chain.results.dc`. Additionally, it
    adjusts the `v_mp` and `p_mp` values by considering `v_pv_wire_drop` and
    `p_pv_wire_loss` from `model_chain.results.dc_ohmic_losses`.

    Returns:
        ModelChain: The modified ModelChain instance with updated `results.ac` based
        on MPPT calculations. Returns `None` if the necessary DC model outputs
        are not present.

    Side Effects:
        - Modifies `model_chain.results.ac` with the calculated MPPT output power.
        - Modifies `model_chain.results.dc_ohmic_losses` when MPPT output is zero.
        - This function does not directly modify input parameters but updates the
          `model_chain` instance's properties based on calculations.

    Notes:
        The results are stored as a pandas Series or a tuple of Series with the same
        index as `model_chain.results.dc`, facilitating integration with further
        PV system analysis or performance calculations.
    """

    #helper function
    def calculate_mppt_output_current(p_mppt_out, v_batt, r_batt_wire):
        if r_batt_wire == 0:
            return p_mppt_out/v_batt
        else:
            # solve the quadratic equation for the current
            return (-1*v_batt + np.sqrt(v_batt**2 + 4*p_mppt_out*r_batt_wire))/(2*r_batt_wire)
    
    def calculate_power(v_mp, v_oc, p_mp):
        p_ac = np.empty_like(p_mp)
        v_batt = model_chain.system.inverter_parameters.get('v_batt', 14.2)
        v_start_delta = v_batt + model_chain.system.inverter_parameters.get('v_start_delta', 5.0)
        v_continue_delta = v_batt + model_chain.system.inverter_parameters.get('v_continue_delta', 1.0)
        v_start = v_batt + v_start_delta
        v_continue = v_batt + v_continue_delta
        mppt_eff = model_chain.system.inverter_parameters.get('mppt_eff', 1.0)
        r_batt_wire = model_chain.system.inverter_parameters.get('r_batt_wire', 0)

        voltage_drops = getattr(model_chain.results.dc_ohmic_losses, 'v_pv_wire_drop', np.zeros_like(v_mp))
        power_losses = getattr(model_chain.results.dc_ohmic_losses, 'p_pv_wire_loss', np.zeros_like(v_mp))

        on = False
        for i, (vmp, voc, pmp, vdrop_in, ploss_in) in enumerate(zip(v_mp, v_oc, p_mp, voltage_drops, power_losses)):
            p_mppt_in = pmp - ploss_in
            p_mppt_out = mppt_eff * p_mppt_in
            i_out = calculate_mppt_output_current(p_mppt_out, v_batt, r_batt_wire)
            ohmic_loss_out = i_out**2 * r_batt_wire
            voltage_drop_out = i_out * r_batt_wire
            v_mppt_in = vmp - vdrop_in
                
            if on:
                if v_mppt_in < v_continue + voltage_drop_out:
                    on = False
                    p_ac[i] = 0
                else:
                    p_ac[i] = p_mppt_out - ohmic_loss_out
            else:
                if voc >= v_start and v_mppt_in >= v_continue + voltage_drop_out:
                    on = True
                    p_ac[i] = p_mppt_out - ohmic_loss_out
                else:
                    p_ac[i] = 0
                
        if isinstance(p_mp, pd.Series):
            p_ac = pd.Series(p_ac, index=p_mp.index)
            p_ac.name = 'p_ac'
        return p_ac
    
    def zero_dc_ohmic_losses_when_p_ac_is_zero(p_ac, dc_ohmic_losses):
        mask = p_ac == 0
        dc_ohmic_losses.loc[mask] = 0
        return dc_ohmic_losses

    # test for required values
    if not all(key in model_chain.results.dc for key in ['v_mp', 'v_oc', 'p_mp']):
        model_chain.results.ac = None
        return model_chain

    if isinstance(model_chain.results.dc, tuple):
        model_chain.results.ac = tuple(calculate_power(df['v_mp'], 
                                                       df['v_oc'], 
                                                       df['p_mp']) for df in model_chain.results.dc)
    else:
        model_chain.results.ac = calculate_power(model_chain.results.dc['v_mp'], 
                                                 model_chain.results.dc['v_oc'], 
                                                 model_chain.results.dc['p_mp'])
    
    # zero out dc ohmic losses when p_ac is zero
    if model_chain.results.dc_ohmic_losses is not None and \
    all(key in model_chain.results.dc_ohmic_losses for key in ['p_pv_wire_loss', 'v_pv_wire_drop']):
        model_chain.results.dc_ohmic_losses = zero_dc_ohmic_losses_when_p_ac_is_zero(
            model_chain.results.ac, model_chain.results.dc_ohmic_losses)
    
    return model_chain

def pv_wire_loss(model_chain):
    """
    Calculate the PV side wire ohmic loss and PV side voltage drop.

    Ohmic power loss is calculated as: 
    model_chain.result.dc['i_mp']**2 * model_chain.system.inverter_parameters['r_pv_wire']

    Ohmic voltage drop is calulated as: 
    model_chain.result.dc['i_mp'] * model_chain.system.inverter_parameters['r_pv_wire']
    
    Parameters:
        model_chain (ModelChain): An instance of a ModelChain that includes the
        PV system configuration and solar position. It must contain the DC model output
        and inverter parameters.

    Optional Parameters:
        - `r_pv_wire` (float): Resistance of PV side wiring, default 0 if not provided.

    The function expects the use of DC models (e.g., "cec" or "desoto") that provide `i_mp` 
    values in `model_chain.results.dc`. 

    Results are stored in model_chain.results.dc_ohmic_loss.

    Returns:
        ModelChain: The modified ModelChain instance with updated `results._dc_ohmic_losses` based
        on the above calculations or `None` if the necessary DC model outputs
        are not present.

    Side Effects:
        - Modifies `model_chain.results.dc_ohmic_losses`
        
    Notes:
        The results are stored as a pandas Dataframes or a tuple of Dataframes with the same
        index as `model_chain.results.dc`, facilitating integration with further
        PV system analysis or performance calculations.
    """
    #TODO We should think about this implementation a bit more.  It assumes ohmic loss when
    # at times the MPPT will not be operating and thus there is no real loss.
    # For our application it might make sense to roll this into the MPPT model.
    # Alternatively we can modify mc.results.dc_ohmic_loss by the MPPT model
    
    #helper function
    def calculate_power_and_voltage_drop(i_mp, r_pv_wire, index):
        p =  i_mp**2 * r_pv_wire
        v = i_mp * r_pv_wire
        points = p, v
        columns = ['p_pv_wire_loss', 'v_pv_wire_drop']
        return pd.DataFrame(np.array(points).T, index=index, columns=columns)
    
    # test for required values
    
    if 'i_mp' not in model_chain.results.dc:
        model_chain.results.dc_ohmic_losses = None
        return model_chain
    
    r_pv_wire = model_chain.system.inverter_parameters.get('r_pv_wire', 0)

    if isinstance(model_chain.results.dc, tuple):
        model_chain.results.dc_ohmic_losses = tuple(calculate_power_and_voltage_drop(
            df['i_mp'], 
            r_pv_wire,
            df['i_mp'].index
            ) for df in model_chain.results.dc)
    else:
        model_chain.results.dc_ohmic_losses = calculate_power_and_voltage_drop(
            model_chain.results.dc['i_mp'], 
             r_pv_wire,
             model_chain.results.dc['i_mp'].index
             )

    return model_chain
    
def getWireResistance(awg, length):
    """
    Calculates the wire resistance of a wire of guage `awg` and length `length`.  Uses
    resistance per 1000ft values from Ancor technical data.

    Parameters
    ----------
    awg : str
        The American Wire Guage (AWG) of the wire. Must be a string of the form "XXAWG" where
        XX is the guage number 4/0 through 18.
    length : float
        The length of the wire in feet.  Include the total length of the wire,
        including the positive and negative conductors.
    """
    awg_r_per_1000ft = {
        "4/0AWG": 0.05,
        "3/0AWG": 0.06,
        "2/0AWG": 0.08,
        "1/0AWG": 0.1,
        "1AWG": 0.13,
        "2AWG": 0.16,
        "4AWG": 0.24,
        "6AWG": 0.4,
        "8AWG": 0.62,
        "10AWG": 1.00,
        "12AWG": 1.75,
        "14AWG": 2.53,
        "16AWG": 4.00,
        "18AWG": 6.48
    }

    if awg not in awg_r_per_1000ft:
        raise ValueError("Invalid wire guage.  Must be one of: " + ", ".join(awg_r_per_1000ft.keys()))  
    return awg_r_per_1000ft[awg] * length / 1000