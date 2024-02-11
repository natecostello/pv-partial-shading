import pvlib
import pvmismatch
from pvmismatch.contrib.gen_coeffs import gen_two_diode
import numpy as np
import matplotlib.pyplot as plt

# TODO figure out if alpha_sc units are correct or inconsistent between functions/modules
# Helper functions for two-diode model parameters

def compare_two_diode_solution_x(module_params, x):
    pvc = pvmismatch.pvcell.PVcell(Isat1_T0=x[0], Isat2_T0=x[1], Rs=x[2], Rsh=x[3], Isc0_T0=module_params['I_sc_ref'], alpha_Isc=module_params['alpha_sc_percent'])
    mpp = np.argmax(pvc.Pcell)  # find the index of the max power point
    specff = module_params['I_mp_ref']*module_params['V_mp_ref'] / module_params['I_sc_ref'] / module_params['V_oc_ref']
    ff = (pvc.Icell[mpp] * pvc.Vcell[mpp]) / pvc.Isc / pvc.Voc
    acc_i_sc = (pvc.Isc - module_params['I_sc_ref']) / module_params['I_sc_ref']
    acc_v_oc = (pvc.Voc * module_params['N_s'] - module_params['V_oc_ref']) / module_params['V_oc_ref']
    acc_i_mp = (pvc.Icell[mpp][0] - module_params['I_mp_ref']) / module_params['I_mp_ref']
    acc_v_mp = (pvc.Vcell[mpp][0]*module_params['N_s'] - module_params['V_mp_ref']) / module_params['V_mp_ref']
    acc_ff = (ff[0] - specff) / specff
    report_string =  '2    Diode: Isc: {:.4}, Voc: {:.4}, Imp: {:.4}, Vmp: {:.4}, FF: {:.4} \n'.format(pvc.Isc, pvc.Voc * module_params['N_s'], pvc.Icell[mpp][0], pvc.Vcell[mpp][0]*module_params['N_s'], ff[0])
    report_string += 'Spec sheet: Isc: {:.4}, Voc: {:.4}, Imp: {:.4}, Vmp: {:.4}, FF: {:.4} \n'.format(module_params['I_sc_ref'], module_params['V_oc_ref'], module_params['I_mp_ref'], module_params['V_mp_ref'], specff)
    report_string += 'Difference: Isc: {:.2%}, Voc: {:.2%}, Imp: {:.2%}, Vmp: {:.2%}, FF: {:.2%}'.format(acc_i_sc, acc_v_oc, acc_i_mp, acc_v_mp, acc_ff)
    return report_string

def compare_two_diode_solution(module_params):
    pvc = pvmismatch.pvcell.PVcell(Isat1_T0=module_params['Isat1_T0'], Isat2_T0=module_params['Isat2_T0'], Rs=module_params['Rs_2d'], Rsh=module_params['Rsh_2d'], Isc0_T0=module_params['I_sc_ref'], alpha_Isc=module_params['alpha_sc_percent'])
    mpp = np.argmax(pvc.Pcell)  # find the index of the max power point
    specff = module_params['I_mp_ref']*module_params['V_mp_ref'] / module_params['I_sc_ref'] / module_params['V_oc_ref']
    ff = (pvc.Icell[mpp] * pvc.Vcell[mpp]) / pvc.Isc / pvc.Voc
    acc_i_sc = (pvc.Isc - module_params['I_sc_ref']) / module_params['I_sc_ref']
    acc_v_oc = (pvc.Voc * module_params['N_s'] - module_params['V_oc_ref']) / module_params['V_oc_ref']
    acc_i_mp = (pvc.Icell[mpp][0] - module_params['I_mp_ref']) / module_params['I_mp_ref']
    acc_v_mp = (pvc.Vcell[mpp][0]*module_params['N_s'] - module_params['V_mp_ref']) / module_params['V_mp_ref']
    acc_ff = (ff[0] - specff) / specff
    report_string = '2 Diode parameters: Isat1_T0: {:.4}, Isat2_T0: {:.4}, Rs: {:.4}, Rsh: {:.4} \n'.format(module_params['Isat1_T0'], module_params['Isat2_T0'], module_params['Rs_2d'], module_params['Rsh_2d'])
    report_string +=  '    2 Diode: Isc: {:.4}, Voc: {:.4}, Imp: {:.4}, Vmp: {:.4}, FF: {:.4} \n'.format(pvc.Isc, pvc.Voc * module_params['N_s'], pvc.Icell[mpp][0], pvc.Vcell[mpp][0]*module_params['N_s'], ff[0])
    report_string += '  Spec sheet: Isc: {:.4}, Voc: {:.4}, Imp: {:.4}, Vmp: {:.4}, FF: {:.4} \n'.format(module_params['I_sc_ref'], module_params['V_oc_ref'], module_params['I_mp_ref'], module_params['V_mp_ref'], specff)
    report_string += '  Difference: Isc: {:.2%}, Voc: {:.2%}, Imp: {:.2%}, Vmp: {:.2%}, FF: {:.2%}'.format(acc_i_sc, acc_v_oc, acc_i_mp, acc_v_mp, acc_ff)
    return report_string

def plot_compare_two_diode_solution_STC(module_params):
    pvc = pvmismatch.pvcell.PVcell(
        Isat1_T0=module_params['Isat1_T0'], 
        Isat2_T0=module_params['Isat2_T0'], 
        Rs=module_params['Rs_2d'], 
        Rsh=module_params['Rsh_2d'], 
        Isc0_T0=module_params['I_sc_ref'], 
        alpha_Isc=module_params['alpha_sc_percent'])
    
    # get index max power point
    mpp = np.argmax(pvc.Pcell)

    # use pvlib to get the full IV curve using CEC model
    params1stc = pvlib.pvsystem.calcparams_cec(
        effective_irradiance=1000,
        temp_cell=module_params['temp_ref'], 
        alpha_sc=module_params['alpha_sc'], 
        a_ref=module_params['a_ref'],
        I_L_ref=module_params['I_L_ref'], 
        I_o_ref=module_params['I_o_ref'], 
        R_sh_ref=module_params['R_sh_ref'],
        R_s=module_params['R_s'], 
        Adjust=module_params['Adjust'])
    
    iv_params1stc = pvlib.pvsystem.singlediode(*params1stc, ivcurve_pnts=100, method='newton')

    # use pvmm to get full IV curve using 2-diode model parameters
    #panel_layout = pvmismatch.pvmodule.standard_cellpos_pat(9, [2, 2])
    panel_layout = module_params['cell_layout']
    
    pvm = pvmismatch.pvmodule.PVmodule(
        cell_pos=panel_layout, 
        Vbypass=[-0.5, -0.5], 
        pvcells=[pvc]*module_params['N_s'])
    

    # make some comparison plots
    pvm.plotMod()  # plot the pvmm module
    plt.tight_layout()

    # get axes for IV curve
    f, ax = plt.gcf(), plt.gca()
    ax0 = f.axes[0]
    ax0.plot(iv_params1stc['v'], iv_params1stc['i'], '--')
    ax0.plot(0, iv_params1stc['i_sc'], 'o', mfc='none', mec='orange')
    ax0.plot(iv_params1stc['v_oc'], 0, 'o', mfc='none', mec='orange')
    ax0.plot(iv_params1stc['v_mp'], iv_params1stc['i_mp'], 'o', mfc='none', mec='orange')
    ax0.set_ylim([0, np.ceil(iv_params1stc['i_sc'].max())])
    ax0.set_xlim([0, np.ceil(iv_params1stc['v_oc'].max())])

    ax0.plot(0, pvm.Isc.mean(), 'o', mfc='none', mec='b')
    ax0.plot(pvm.Voc.sum(), 0, 'o', mfc='none', mec='b')

    mpp = np.argmax(pvm.Pmod)
    #ax0.plot(pvm.Vcell[mpp], mpp.Icell[mpp], '_k')
    ax0.plot(pvm.Vmod[mpp], pvm.Imod[mpp], 'o', mfc='none', mec='b')

    iv_params1stc['p'] = iv_params1stc['v'] * iv_params1stc['i']
    ax1 = f.axes[1]
    ax1.plot(iv_params1stc['v'], iv_params1stc['p'], '--')

    ax1.plot(iv_params1stc['v_mp'], iv_params1stc['p_mp'], 'o', mfc='none', mec='orange')
    ax1.plot(pvm.Vmod[mpp], pvm.Pmod[mpp], 'o', mfc='none', mec='b')
    ax1.set_xlim([0, np.ceil(iv_params1stc['v_oc'].max())])

def plot_compare_two_diode_solution(module_params, effective_irradiance, temp_cell):
    pvc = pvmismatch.pvcell.PVcell(
        Isat1_T0=module_params['Isat1_T0'], 
        Isat2_T0=module_params['Isat2_T0'], 
        Rs=module_params['Rs_2d'], 
        Rsh=module_params['Rsh_2d'], 
        Isc0_T0=module_params['I_sc_ref'], 
        alpha_Isc=module_params['alpha_sc_percent'])
    
    # get index max power point
    mpp = np.argmax(pvc.Pcell)

    # use pvlib to get the full IV curve using CEC model
    params1stc = pvlib.pvsystem.calcparams_cec(
        effective_irradiance=effective_irradiance,
        temp_cell=temp_cell, 
        alpha_sc=module_params['alpha_sc'], 
        a_ref=module_params['a_ref'],
        I_L_ref=module_params['I_L_ref'], 
        I_o_ref=module_params['I_o_ref'], 
        R_sh_ref=module_params['R_sh_ref'],
        R_s=module_params['R_s'], 
        Adjust=module_params['Adjust'])
    
    iv_params1stc = pvlib.pvsystem.singlediode(*params1stc, ivcurve_pnts=100, method='newton')

    # use pvmm to get full IV curve using 2-diode model parameters
    #panel_layout = pvmismatch.pvmodule.standard_cellpos_pat(9, [2, 2])
    panel_layout = module_params['cell_layout']
    
    pvm = pvmismatch.pvmodule.PVmodule(
        cell_pos=panel_layout, 
        Vbypass=[-0.5, -0.5], 
        pvcells=[pvc]*module_params['N_s'])
    
    pvm.setSuns(effective_irradiance/1000)
    pvm.setTemps(temp_cell+273.15)

    # make some comparison plots
    pvm.plotMod()  # plot the pvmm module
    plt.tight_layout()

    # get axes for IV curve
    f, ax = plt.gcf(), plt.gca()
    ax0 = f.axes[0]
    ax0.plot(iv_params1stc['v'], iv_params1stc['i'], '--')
    ax0.plot(0, iv_params1stc['i_sc'], 'o', mfc='none', mec='orange')
    ax0.plot(iv_params1stc['v_oc'], 0, 'o', mfc='none', mec='orange')
    ax0.plot(iv_params1stc['v_mp'], iv_params1stc['i_mp'], 'o', mfc='none', mec='orange')
    ax0.set_ylim([0, np.ceil(iv_params1stc['i_sc'].max())])
    ax0.set_xlim([0, np.ceil(iv_params1stc['v_oc'].max())])

    ax0.plot(0, pvm.Isc.mean(), 'o', mfc='none', mec='b')
    ax0.plot(pvm.Voc.sum(), 0, 'o', mfc='none', mec='b')

    mpp = np.argmax(pvm.Pmod)
    #ax0.plot(pvm.Vcell[mpp], mpp.Icell[mpp], '_k')
    ax0.plot(pvm.Vmod[mpp], pvm.Imod[mpp], 'o', mfc='none', mec='b')

    iv_params1stc['p'] = iv_params1stc['v'] * iv_params1stc['i']
    ax1 = f.axes[1]
    ax1.plot(iv_params1stc['v'], iv_params1stc['p'], '--')

    ax1.plot(iv_params1stc['v_mp'], iv_params1stc['p_mp'], 'o', mfc='none', mec='orange')
    ax1.plot(pvm.Vmod[mpp], pvm.Pmod[mpp], 'o', mfc='none', mec='b')
    ax1.set_xlim([0, np.ceil(iv_params1stc['v_oc'].max())])


def calc_two_diode_params(i_sc_ref, v_oc_ref, i_mp_ref, v_mp_ref, n_s, parallel_strings, temp_ref=25):
    """
    Calculate the two-diode model parameters from module specification parameters.

    Args:
        i_sc_ref (float): short-circuit current at reference conditions, A
        v_oc_ref (float): open-circuit voltage at reference conditions, V
        i_mp_ref (float): current at maximum power point at reference conditions, A
        v_mp_ref (float): voltage at maximum power point at reference conditions, V
        n_s (int): number of cells in series
        parallel_strings (int): number of parallel strings
        temp_ref (float): reference temperature, C

    Returns:
        tuple: (Isat1_T0, Isat2_T0, Rs, Rsh)
    """
    def last_guess(sol):
        isat1 = np.exp(sol.x[0])
        isat2 = np.exp(sol.x[1])
        rs = sol.x[2] ** 2.0
        rsh = sol.x[3] ** 2.0
        return isat1, isat2, rs, rsh

    args = (i_sc_ref, v_oc_ref, i_mp_ref, v_mp_ref, n_s, parallel_strings, temp_ref)

    # First Calculation
    x, sol = gen_two_diode(*args)

    # Iterate
    x = last_guess(sol)
    x, sol = gen_two_diode(*args, x0=x)
    x = last_guess(sol)
    x, sol = gen_two_diode(*args, x0=x)
    x = last_guess(sol)
    x, sol = gen_two_diode(*args, x0=x)
    x = last_guess(sol)
    x, sol = gen_two_diode(*args, x0=x)

    return x

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
        'alpha_sc_percent': 0.0009, # 1/C = 0.0005*1.8
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
        'Adjust': None,
        'Parallel_strings': 1
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

    Isat1_T0, Isat2_T0, Rs, Rsh = calc_two_diode_params(
        module_parameters['I_sc_ref'],
        module_parameters['V_oc_ref'],
        module_parameters['I_mp_ref'],
        module_parameters['V_mp_ref'],
        module_parameters['N_s'],
        module_parameters['Parallel_strings'],
        module_parameters['temp_ref']
        )

    module_parameters['Isat1_T0'] = Isat1_T0
    module_parameters['Isat2_T0'] = Isat2_T0
    module_parameters['Rs_2d'] = Rs
    module_parameters['Rsh_2d'] = Rsh
    module_parameters['cell_layout'] = pvmismatch.pvmodule.standard_cellpos_pat(18, [2, 2])


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
        'alpha_sc_percent': 0.0009, # 1/C = 0.0005*1.8  
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
        'Adjust': None,
        'Parallel_strings': 1
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

    Isat1_T0, Isat2_T0, Rs, Rsh = calc_two_diode_params(
        module_parameters['I_sc_ref'],
        module_parameters['V_oc_ref'],
        module_parameters['I_mp_ref'],
        module_parameters['V_mp_ref'],
        module_parameters['N_s'],
        module_parameters['Parallel_strings'],
        module_parameters['temp_ref']
        )

    module_parameters['Isat1_T0'] = Isat1_T0
    module_parameters['Isat2_T0'] = Isat2_T0
    module_parameters['Rs_2d'] = Rs
    module_parameters['Rsh_2d'] = Rsh
    module_parameters['cell_layout'] = pvmismatch.pvmodule.standard_cellpos_pat(9, [2, 2])

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
        'N_s': 34,  # 2 strings of 34 cells 
        'I_sc_ref': 11.05, 
        'V_oc_ref': 23.0, 
        'I_mp_ref': 10.42, 
        'V_mp_ref': 19.2, 
        'alpha_sc': 0.009945, # A/C = 0.0005*11.05*1.8
        'alpha_sc_percent': 0.0009, # 1/C = 0.0005*1.8
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
        'Adjust': None,
        'Parallel_strings': 2
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

    Isat1_T0, Isat2_T0, Rs, Rsh = calc_two_diode_params(
        module_parameters['I_sc_ref'],
        module_parameters['V_oc_ref'],
        module_parameters['I_mp_ref'],
        module_parameters['V_mp_ref'],
        module_parameters['N_s'],
        module_parameters['Parallel_strings']
        )

    module_parameters['Isat1_T0'] = Isat1_T0
    module_parameters['Isat2_T0'] = Isat2_T0
    module_parameters['Rs_2d'] = Rs
    module_parameters['Rsh_2d'] = Rsh
    module_parameters['cell_layout'] = pvmismatch.pvmodule.standard_cellpos_pat(17, [1, 1])

    return module_parameters