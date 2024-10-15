# class for optimization functions which is called in element_tracker.py

import numpy as np
from scipy import constants
from scipy import stats

class Optimisation_Functions:
    def __init__(self, tracked_beam, settings):
        self.tracked_beam = tracked_beam
        self.nbins = 100
        self.settings = settings

    # switch case for deciding the optimisation method + passing the optimisation function back
    def get_opt_func_val(self):
        if isinstance(self.settings, list):
            self.nopt = len(self.settings)
            # write here method for getting the multi-objective run working
            # one-liner for the return below can be used in the function for getting multi-runs
            # make a function to add in quadrature all of the opt functions
            # set position
            #! need to modify the optimisation functions for a self.settings[i] case! for penalty values 
            obj_func_spec = [[],[],[]]
            for i in range(self.nopt):
                try:
                    obj_func_spec[0].append(self.settings[i].get('obj_func'))
                    if self.settings[i].get('pos') != None:
                        obj_func_spec[1].append(self.settings[i].get('pos'))
                    else:
                        obj_func_spec[1].append(-1)
                    obj_func_spec[2].append(self.settings[i].get('p_val'))
                except AttributeError:
                    print("Enter correct optimisation functions")
            # return the added in quadrature penalty value
            return np.sqrt(sum([getattr(self, obj_func_spec[0][i] + "_opt")(obj_func_spec[1][i], obj_func_spec[2][i]) for i in range(self.nopt)]))
        else:
            # set -1 i.e. end of beamline is an optional position isn't passed
            if self.settings.get('pos') != None:
                pos = self.settings.get('pos')
            else:
                pos = -1
            # selection + running of optmisation function
            try:
                return getattr(self, self.settings.get('obj_func') + "_opt")(pos, self.settings.get('p_val'))
            except AttributeError:
                print("Enter a correct optimisation function")
            
    # optimizing for minimum bunch duration or to specific value at FEL
    def sig_t_opt(self, pos, p_val=None):
        if p_val != None:
            return abs(p_val/constants.femto - np.sqrt(np.mean(self.tracked_beam[pos][0]**2))/(constants.c*constants.femto))
        else:
            return np.sqrt(np.mean(self.tracked_beam[pos][0]**2))/(constants.c*constants.femto)
        
    # optimizing for minimum full width bunch duration or to specific value at FEL
    def fw_t_opt(self, pos, p_val=None):
        if p_val != None:
            return abs(p_val/constants.femto - (max(self.tracked_beam[pos][0])-min(self.tracked_beam[pos][0]))/(constants.c*constants.femto)) 
        else:
            return (max(self.tracked_beam[pos][0])-min(self.tracked_beam[pos][0]))/(constants.c*constants.femto)

    # optimizing for minimum energy spread or to specific value at FEL
    def sig_e_opt(self, pos, p_val=None):
        if p_val != None:
            return abs(p_val - np.sqrt(np.mean(self.tracked_beam[pos][1]**2)))
        else:
            return np.sqrt(np.mean(self.tracked_beam[pos][1]**2))
        
    # optimizing for maximum peak current
    def Ipk_opt(self, pos, p_val=None):
        counts, bins = np.histogram(self.tracked_beam[pos][0]/(constants.c*constants.femto), self.nbins)
        Ipk = max(counts)*(self.Q*constants.pico*constants.c*self.nbins)/(sum(counts)*(max(self.tracked_beam[pos][0])-min(self.tracked_beam[pos][0])))
        if p_val != None:
            return abs(p_val - Ipk)
        else:
            # minimization algorithm ...
            return -Ipk

    # optimising for a linear phase space
    #! standard optimisation is making linear as possible
    #! penalty function is pearson_correlation > p_val     
    def lin_opt(self, pos, p_val=None):
        if p_val != None:
            return (1 - abs(p_val)) - (1 - abs(stats.pearsonr(self.tracked_beam[pos][0], self.tracked_beam[pos][1], alternative = 'two-sided').statistic))
        else:
            return 1 - abs(stats.pearsonr(self.tracked_beam[pos][0], self.tracked_beam[pos][1], alternative = 'two-sided').statistic)