# class for optimization functions which is called in element_tracker.py

import numpy as np
from scipy import constants

class Optimisation_Functions:
    def __init__(self, tracked_beam, settings):
        self.tracked_beam = tracked_beam
        self.nbins = 100
        self.settings = settings
        # set -1 i.e. end of beamline is an optional position isn't passed
        if self.settings.get('pos') != None:
            self.pos = self.settings.get('pos')
        else:
            self.pos = -1

    # optimizing for minimum bunch duration or to specific value at FEL
    def sig_t_opt(self):
        if self.settings.get('p_vals') != None:
            if self.settings.get('p_vals').get('p_sigt') != None:
                return abs(self.settings.get('p_vals').get('p_sigt')/constants.femto - np.sqrt(np.mean(self.tracked_beam[self.pos][0]**2))/(constants.c*constants.femto))
            else:
               return np.sqrt(np.mean(self.tracked_beam[self.pos][0]**2))/(constants.c*constants.femto) 
        else:
            return np.sqrt(np.mean(self.tracked_beam[self.pos][0]**2))/(constants.c*constants.femto)

    # optimizing for minimum energy spread or to specific value at FEL
    def sig_e_opt(self):
        if self.settings.get('p_vals') != None:
            if self.settings.get('p_vals').get('p_sige') != None:
                return abs(self.settings.get('p_vals').get('p_sige') - np.sqrt(np.mean(self.tracked_beam[self.pos][1]**2)))
            else:
                return np.sqrt(np.mean(self.tracked_beam[self.pos][1]**2))
        else:
            return np.sqrt(np.mean(self.tracked_beam[self.pos][1]**2))
        
        # optimizing for maximum peak current
    def Ipk_opt(self, Q):
        counts, bins = np.histogram(self.tracked_beam[self.pos][0]/(constants.c*constants.femto), self.nbins)
        Ipk = -max(counts)*(Q*constants.pico*constants.c*self.nbins)/(sum(counts)*(max(self.tracked_beam[self.pos][0])-min(self.tracked_beam[self.pos][0])))
        print(Ipk)
        return Ipk
        
    def sig_t_e_opt(self):
        return np.sqrt(self.sig_e_opt(self.pos)**2 + self.sig_t_opt(self.pos)**2)