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

#---------------------
# Beamline Definition
#---------------------

SHINE_init = beam_generation.Beam_Generation("param", 100, 1e4, 6, 0.8638, 8.6e-4, 100.6e6)
SHINE_beam = SHINE_init.gaussian_bunch()

SHINE_dict = munch.munchify({'L1': {'type': 'RF', 'f': 1300e6, 'V': 231e6, 'phi': np.radians(12.1)}, 'HC3': {'type': 'RF', 'f': 3900e6, 'V': 68e6, 'phi': np.radians(153)}, 'BC1': {'type': 'BC', 'R56': -62e-3} ,'L2': {'type': 'RF', 'f': 1300e6, 'V': 2168e6, 'phi': np.radians(25)}, 'BC2': {'type': 'BC', 'R56': -43.95e-3}, 'L3': {'type': 'RF', 'f': 1300e6, 'V': 2891e6, 'phi': 0}, 'L4': {'type': 'RF', 'f': 1300e6, 'V': 3614e6, 'phi': 0}})

SHINE_tracker = element_tracker.Element_Track(SHINE_beam, 100, SHINE_dict)
SHINE_tracker.synchrotron_radiation(20) 
#--------------
# Standard Run
#--------------

#SHINE_track = SHINE_tracker.track_1D()

#SHINE_plt = comp_plt.Track_Plot(SHINE_dict, SHINE_beam, SHINE_track, "SHINE", SHINE_init.bunch_charge, SHINE_init.npart)
#E_mean = SHINE_plt.beam_energy_plt()
#sig_t = SHINE_plt.beam_sigt_plt()
#sig_e = SHINE_plt.beam_sige_plt()
#Ipk = SHINE_plt.current_plt()
#print(Ipk)
#--------------
# Optimisation
#--------------

# the optimization can be limited to certain variables by specifying the bounds as identical to the values 
#SHINE_bounds = [(1300e6, 1300e6), (231e6, 231e6), (-np.pi, np.pi), (3900e6, 3900e6), (68e6, 68e6), (-np.pi, np.pi), (-10e-2, 10e-2), (1300e6, 1300e6), (2168e6, 2168e6), (np.radians(25), np.radians(25)), (-43.95e-3, -43.95e-3), (1300e6, 1300e6), (2891e6, 2891e6), (0, 0), (1300e6, 1300e6), (3614e6, 3614e6), (0, 0)]
#SHINE_bounds = [(1300e6, 1300e6), (231e6, 231e6), (-np.pi, np.pi), (3900e6, 3900e6), (68e6, 68e6), (-np.pi, np.pi), (-10e-2, 10e-2), (1300e6, 1300e6), (2168e6, 2168e6), (-np.pi, np.pi), (-10e-2, 10e-2), (1300e6, 1300e6), (2891e6, 2891e6), (-np.pi, np.pi), (1300e6, 1300e6), (3614e6, 3614e6), (-np.pi, np.pi)]

# penalty function values can be passed as a dictionary
#pvals = {'p_sigt': 20e-15, 'p_sige': 1e-5}

# selected objective function
#obj_func = 'sig_t'
#obj_func = 'Ipk'

# if no penalty function values are passed in the arguments then minimises the objective function
#result = minimize(SHINE_tracker.opt_func, SHINE_tracker.initial_guess(), args = (obj_func), method = 'Nelder-Mead', bounds = SHINE_bounds, options = {'maxiter':2e4, 'fatol':0.01})

#print('----------------')
#print('Optimisation Results:')
#print("penalty function/result: ", result.fun)
#SHINE_init_param = np.array([list(SHINE_bounds[i]) for i in range(len(SHINE_bounds))])[:,0]

#param_names = []
#for k,v in SHINE_dict.items():
#    param_names.append([k + '_' + s for s in list(v.keys())])
#param_names = np.array([i for i in chain.from_iterable(param_names) if 'type' not in i])

#SHINE_opt_param = np.array(result.x)
#index_opt_param = np.flatnonzero((~np.equal(SHINE_init_param, SHINE_opt_param)).astype(int))
#val_opt_param = SHINE_opt_param[index_opt_param]
#name_opt_param = param_names[index_opt_param] 
#print("Optimised Variables")
#for i in range(len(name_opt_param)):
#    print(name_opt_param[i] + ": " + str(val_opt_param[i])) 
#if result.success == False:
#    print(result.message)
#print('----------------')

# run the solution
#SHINE_track_opt = SHINE_tracker.track_1D()

#-------
# PLOTS
#-------
#SHINE_plt = comp_plt.Track_Plot(SHINE_dict, SHINE_beam, SHINE_track_opt, "SHINE", SHINE_init.bunch_charge, SHINE_init.npart, save_plots=True)

#SHINE_plt.beam_sigt_sige_phase_space()
#Ipk_final = SHINE_plt.current_plt()

#E_mean = SHINE_plt.beam_energy_plt()
#sig_t = SHINE_plt.beam_sigt_plt()
#sig_e = SHINE_plt.beam_sige_plt()
#Ipk = SHINE_plt.beam_Ipk_plt()
#print(sig_t)
#print(Ipk)
#plt.show()


