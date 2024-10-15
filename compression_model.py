# compression 1D model

import numpy as np
import munch
from scipy import constants
from matplotlib import pyplot as plt
import beam_generation
import element_tracker
import compression_plotting as comp_plt

from itertools import chain
from scipy.optimize import minimize

if __name__ == "__main__":

    #---------------------
    # Beamline Definition
    #---------------------

    #SHINE_init = beam_generation.Beam_Generation("param", 100, 1e4, 6, 0.8638, 8.6e-4, 100.6e6)
    #SHINE_beam = SHINE_init.gaussian_bunch()

    #SHINE_dict = munch.munchify({'L1': {'type': 'RF', 'f': 1300e6, 'V': 231e6, 'phi': np.radians(12.1)}, 'HC3': {'type': 'RF', 'f': 3900e6, 'V': 68e6, 'phi': np.radians(153)}, 'BC1': {'type': 'BC', 'R56': -62e-3} ,'L2': {'type': 'RF', 'f': 1300e6, 'V': 2168e6, 'phi': np.radians(25)}, 'BC2': {'type': 'BC', 'R56': -43.95e-3}, 'L3': {'type': 'RF', 'f': 1300e6, 'V': 2891e6, 'phi': 0}, 'L4': {'type': 'RF', 'f': 1300e6, 'V': 3614e6, 'phi': 0}})

    #SHINE_tracker = element_tracker.Element_Track(SHINE_beam, 100, SHINE_dict)
    #SHINE_tracker.synchrotron_radiation(20) 

    # parameter based generation
    #UKXFEL_init = beam_generation.Beam_Generation("param", 150, 1e4, 6, 6, 2.5e-4, 80e6)
    #UKXFEL_beam = UKXFEL_init.cubic_sige_bunch()

    # from file
    UKXFEL_init = beam_generation.Beam_Generation("file","VHF_Linearisation.h5")
    UKXFEL_beam = UKXFEL_init.get_file_LPS_data()

    #UKXFEL_dict = munch.munchify({'CAV1': {'type': 'RF', 'f': 650e6, 'V': 30e6, 'phi': -0.31542948109941077}, 'BC1': {'type': 'BC', 'R56': 0.48195151982955886}, 'HC3': {'type': 'RF', 'f': 3900e6, 'V': 30e6, 'phi': np.radians(-45)}, 'CAV2': {'type': 'RF', 'f': 1300e6, 'V': 200e6, 'phi': np.radians(7)},'BC2': {'type': 'BC', 'R56': -0.09}, 'CAV3': {'type': 'RF', 'f': 1300e6, 'V': 400e6, 'phi': np.radians(-7)}}) 
    UKXFEL_dict = munch.munchify({'CAV1': {'type': 'RF', 'f': 650e6, 'V': 30e6, 'phi': -0.31542948109941077}, 'BC1': {'type': 'BC', 'R56': 0.48195151982955886}, 'HC3': {'type': 'RF', 'f': 3900e6, 'V': 30e6, 'phi': np.radians(-45)}})

    UKXFEL_tracker = element_tracker.Element_Track(UKXFEL_beam, UKXFEL_init.bunch_charge, UKXFEL_dict)

    #--------------
    # Standard Run
    #--------------

    #SHINE_track = SHINE_tracker.track_1D()

    #SHINE_plt = comp_plt.Track_Plot(SHINE_dict, SHINE_beam, SHINE_track, "SHINE", SHINE_init.bunch_charge, SHINE_init.npart)
    #E_mean = SHINE_plt.beam_energy_plt()
    #sig_t = SHINE_plt.beam_sigt_plt()
    #sig_e = SHINE_plt.beam_sige_plt()
    #Ipk = SHINE_plt.current_plt()
    #SHINE_plt.beam_sigt_sige_phase_space()

    #UKXFEL_track = UKXFEL_tracker.track_1D()

    #UKXFEL_plt = comp_plt.Track_Plot(UKXFEL_dict, UKXFEL_beam, UKXFEL_track, "UKXFEL", UKXFEL_init.bunch_charge, UKXFEL_init.npart)
    #E_mean = UKXFEL_plt.beam_energy_plt()
    #sig_t = UKXFEL_plt.beam_sigt_plt()
    #sig_e = UKXFEL_plt.beam_sige_plt()
    #Ipk = UKXFEL_plt.current_plt()
    #for i in range(len(UKXFEL_dict.keys())+1):
    #    UKXFEL_plt.beam_sigt_sige_phase_space(i)
    #print(sig_t)
    #print(E_mean)


    #--------------
    # Optimisation
    #--------------

    # the optimization can be limited to certain variables by specifying the bounds as identical to the values 
    #SHINE_bounds = [(1300e6, 1300e6), (231e6, 231e6), (-np.pi, np.pi), (3900e6, 3900e6), (68e6, 68e6), (-np.pi, np.pi), (-10e-2, 10e-2), (1300e6, 1300e6), (2168e6, 2168e6), (np.radians(25), np.radians(25)), (-43.95e-3, -43.95e-3), (1300e6, 1300e6), (2891e6, 2891e6), (0, 0), (1300e6, 1300e6), (3614e6, 3614e6), (0, 0)]

    #UKXFEL_bounds = [(650e6, 650e6), (30e6, 30e6), (-0.31542948109941077, -0.31542948109941077), (0.48195151982955886, 0.48195151982955886), (3900e6, 3900e6), (30e6, 30e6), (np.radians(-60), np.radians(-20)), (1300e6, 1300e6), (200e6, 200e6), (-np.pi, np.pi), (-1, 1), (1300e6, 1300e6), (400e6, 400e6), (np.radians(-7), np.radians(-7))]
    UKXFEL_bounds = [(650e6, 650e6), (30e6, 30e6), (-np.pi, np.pi), (-1.64, 1.64), (3900e6, 3900e6), (30e6, 30e6), (np.radians(-60), np.radians(-20))]

    # settings values optimisation type, position & penalty function values passed as a dictionary
    # position and penalty vals are optional {'obj_func': <type>, <'pos': <pos>>, <'p_vals':{'p_sig': <pen_val>}>}
    settings = [{'obj_func': 'fw_t', 'p_val': 10e-12}, {'obj_func': 'lin', 'p_val': 0.9}]
    #settings = {'obj_func': 'lin', 'pos': -1}

    # if no penalty function values are passed in the arguments then minimises the objective function
    result = minimize(UKXFEL_tracker.opt_func, UKXFEL_tracker.initial_guess(), args = (settings), method = 'Nelder-Mead', bounds = UKXFEL_bounds, options = {'maxiter':2e4, 'fatol':0.01})

    print('----------------')
    print('Optimisation Results:')
    print("penalty function/result: ", result.fun)
    UKXFEL_init_param = np.array([list(UKXFEL_bounds[i]) for i in range(len(UKXFEL_bounds))])[:,0]

    param_names = []
    for k,v in UKXFEL_dict.items():
        param_names.append([k + '_' + s for s in list(v.keys())])
    param_names = np.array([i for i in chain.from_iterable(param_names) if 'type' not in i])

    UKXFEL_opt_param = np.array(result.x)
    index_opt_param = np.flatnonzero((~np.equal(UKXFEL_init_param, UKXFEL_opt_param)).astype(int))
    val_opt_param = UKXFEL_opt_param[index_opt_param]
    name_opt_param = param_names[index_opt_param] 
    print("Optimised Variables")
    for i in range(len(name_opt_param)):
        print(name_opt_param[i] + ": " + str(val_opt_param[i])) 
    if result.success == False:
        print(result.message)
    print('----------------')

    # run the solution
    UKXFEL_track_opt = UKXFEL_tracker.track_1D()

    #-------
    # PLOTS
    #-------
    UKXFEL_plt = comp_plt.Track_Plot(UKXFEL_dict, UKXFEL_beam, UKXFEL_track_opt, "UKXFEL", UKXFEL_init.bunch_charge, UKXFEL_init.npart, save_plots=True)

    for i in range(len(UKXFEL_dict.keys())+1):
        UKXFEL_plt.beam_sigt_sige_phase_space(i)

    Ipk_final = UKXFEL_plt.current_plt()

    E_mean = UKXFEL_plt.beam_energy_plt()
    sig_t = UKXFEL_plt.beam_sigt_plt()
    sig_e = UKXFEL_plt.beam_sige_plt()
    Ipk = UKXFEL_plt.beam_Ipk_plt()
    linearity = UKXFEL_plt.beam_linearity_plt()
    
    print("Rms Bunch Duration: ", sig_t)
    print("Rms Energy Spread: ", sig_e)
    print("Mean Electron Energy: ", E_mean)
    print("Peak Current: ", Ipk)
    print("Pearson Correlation Coefficient: ", linearity)

    plt.show()


