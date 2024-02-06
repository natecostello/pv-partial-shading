import pvlib

RICH_SOLAR_12V = {
    'Name': 'Rich Solar 12V', 
    'Technology': 'Mono-c-Si', 
    'Bifacial': 0, 
    'STC': 200.0862, 
    'PTC': 178.5, 
    'A_c': 2.554, 
    'Length': 1.58, 
    'Width': 0.808, 
    'N_s': 36.0, 
    'I_sc_ref': 11.32, 
    'V_oc_ref': 22.81, 
    'I_mp_ref': 10.74, 
    'V_mp_ref': 18.63, 
    'alpha_sc': 0.006962, 
    'beta_oc': -0.080428, 
    'gamma_r': -0.4567, 
    'T_NOCT': 47.6, 
    'irrad_ref': 1000, 
    'temp_ref': 25, 
    'a_ref': 0.973804158972974, 
    'I_L_ref': 11.44031204496401, 
    'I_o_ref': 7.608255760344541e-10, 
    'R_s': 0.12270248066631148, 
    'R_sh_ref': 197.25444549468273, 
    'Adjust': 7.669332700388107
    }

def get_rich_solar_MEGA200Max_24V():
    module_parameters = {
        'Name': 'Rich Solar MEGA200Max', 
        'Technology': 'Mono-c-Si', 
        'Bifacial': 0, 
        'STC': 200.03, # calculated
        'PTC': None, 
        'A_c': 2.188,  # calculated (not sure if this is total cell area)
        'Length': 1.49, 
        'Width': 0.681, 
        'N_s': 72, 
        'I_sc_ref': 5.83, 
        'V_oc_ref': 45.4, 
        'I_mp_ref': 5.32, 
        'V_mp_ref': 37.6, 
        'alpha_sc': 0.005247, # A/C = 0.0005*5.83*1.8
        'beta_oc': -0.23699,  # V/C = -0.0029*45.4*1.8
        'gamma_r': -0.00702,  # %/C = -0.0039*1.8
        'T_NOCT': 45, 
        'irrad_ref': 1000, 
        'temp_ref': 25, 
        'a_ref': None, 
        'I_L_ref': None, 
        'I_o_ref': None, 
        'R_s': None, 
        'R_sh_ref': None, 
        'Adjust': None
    }

    # Caculate the remaining parameters
    I_L_ref, I_o_ref, R_s, R_sh_ref, a_ref, Adjust = pvlib.ivtools.sdm.fit_cec_sam(
                                                                            celltype='monoSi',
                                                                            v_mp=module_parameters['V_mp_ref'],
                                                                            i_mp=module_parameters['I_mp_ref'],
                                                                            v_oc=module_parameters['V_oc_ref'],
                                                                            i_sc=module_parameters['I_sc_ref'],
                                                                            alpha_sc=module_parameters['alpha_sc'],
                                                                            beta_voc=module_parameters['beta_oc'],
                                                                            gamma_pmp=module_parameters['gamma_r'],
                                                                            cells_in_series=module_parameters['N_s']
                                                                            )
    
    module_parameters['a_ref'] = a_ref
    module_parameters['I_L_ref'] = I_L_ref
    module_parameters['I_o_ref'] = I_o_ref
    module_parameters['R_s'] = R_s
    module_parameters['R_sh_ref'] = R_sh_ref
    module_parameters['Adjust'] = Adjust

    return module_parameters

def get_rich_solar_MEGA200_12V():
    module_parameters = {
        'Name': 'Rich Solar MEGA200', 
        'Technology': 'Mono-c-Si', 
        'Bifacial': 0, 
        'STC': 200.0862, 
        'PTC': 178.5, 
        'A_c': 2.188,  # calculated (not sure if this is total cell area)
        'Length': 1.49, 
        'Width': 0.681, 
        'N_s': 36, 
        'I_sc_ref': 10.2, 
        'V_oc_ref': 24.3, 
        'I_mp_ref': 9.8, 
        'V_mp_ref': 20.4, 
        'alpha_sc': 0.00918,  # A/C = 0.0005*10.2*1.8
        'beta_oc': -0.12685,  # V/C = -0.0029*24.3*1.8
        'gamma_r': -0.00702,  # %/C = -0.0039*1.8
        'T_NOCT': 45, 
        'irrad_ref': 1000, 
        'temp_ref': 25, 
        'a_ref': None, 
        'I_L_ref': None, 
        'I_o_ref': None, 
        'R_s': None, 
        'R_sh_ref': None, 
        'Adjust': None
    }

    # Caculate the remaining parameters
    I_L_ref, I_o_ref, R_s, R_sh_ref, a_ref, Adjust = pvlib.ivtools.sdm.fit_cec_sam(
                                                                            celltype='monoSi',
                                                                            v_mp=module_parameters['V_mp_ref'],
                                                                            i_mp=module_parameters['I_mp_ref'],
                                                                            v_oc=module_parameters['V_oc_ref'],
                                                                            i_sc=module_parameters['I_sc_ref'],
                                                                            alpha_sc=module_parameters['alpha_sc'],
                                                                            beta_voc=module_parameters['beta_oc'],
                                                                            gamma_pmp=module_parameters['gamma_r'],
                                                                            cells_in_series=module_parameters['N_s']
                                                                            )
    
    module_parameters['a_ref'] = a_ref
    module_parameters['I_L_ref'] = I_L_ref
    module_parameters['I_o_ref'] = I_o_ref
    module_parameters['R_s'] = R_s
    module_parameters['R_sh_ref'] = R_sh_ref
    module_parameters['Adjust'] = Adjust

    return module_parameters

def get_renogy_RSP200D_12V():
    module_parameters = {
        'Name': 'Renogy RSP200D', 
        'Technology': 'Mono-c-Si', 
        'Bifacial': 0, 
        'STC': 200.06, 
        'PTC': None, 
        'A_c': 2.188,  # calculated (not sure if this is total cell area)
        'Length': 1.49, 
        'Width': 0.699, 
        'N_s': 68, 
        'I_sc_ref': 11.05, 
        'V_oc_ref': 23.0, 
        'I_mp_ref': 10.42, 
        'V_mp_ref': 19.2, 
        'alpha_sc': 0.009945, # A/C = 0.0005*11.05*1.8
        'beta_oc': -0.11592,  # V/C = -0.0028*23.0*1.8
        'gamma_r': -0.00666,  # %/C = -0.0037*1.8
        'T_NOCT': 47, 
        'irrad_ref': 1000, 
        'temp_ref': 25, 
        'a_ref': None, 
        'I_L_ref': None, 
        'I_o_ref': None, 
        'R_s': None, 
        'R_sh_ref': None, 
        'Adjust': None
    }

    # Caculate the remaining parameters
    I_L_ref, I_o_ref, R_s, R_sh_ref, a_ref, Adjust = pvlib.ivtools.sdm.fit_cec_sam(
                                                                            celltype='monoSi',
                                                                            v_mp=module_parameters['V_mp_ref'],
                                                                            i_mp=module_parameters['I_mp_ref'],
                                                                            v_oc=module_parameters['V_oc_ref'],
                                                                            i_sc=module_parameters['I_sc_ref'],
                                                                            alpha_sc=module_parameters['alpha_sc'],
                                                                            beta_voc=module_parameters['beta_oc'],
                                                                            gamma_pmp=module_parameters['gamma_r'],
                                                                            cells_in_series=module_parameters['N_s']
                                                                            )
    
    module_parameters['a_ref'] = a_ref
    module_parameters['I_L_ref'] = I_L_ref
    module_parameters['I_o_ref'] = I_o_ref
    module_parameters['R_s'] = R_s
    module_parameters['R_sh_ref'] = R_sh_ref
    module_parameters['Adjust'] = Adjust

    return module_parameters