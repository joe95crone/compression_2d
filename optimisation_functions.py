# class for optimization functions which is called in element_tracker.py

import numpy as np
from scipy import constants

class Optimisation_Functions:
    def __init__(self, tracked_beam, *pvals):
        self.tracked_beam = tracked_beam
        self.nbins = 100
        try:
            self.pvals = pvals[0]
        except IndexError:
            pass

    # optimizing for minimum bunch duration or to specific value at FEL
    def sig_t_opt(self):
        try:
            p_sigt = self.pvals.get('p_sigt')
            return abs((p_sigt/constants.femto) - np.sqrt(np.mean(self.tracked_beam[-1][0]**2))/(constants.c*constants.femto))
        except AttributeError:
            return np.sqrt(np.mean(self.tracked_beam[-1][0]**2))/(constants.c*constants.femto)
    
    # optimizing for maximum peak current
    def Ipk_opt(self, Q):
        counts, bins = np.histogram(self.tracked_beam[-1][0]/(constants.c*constants.femto), self.nbins)
        Ipk = -max(counts)*(Q*constants.pico*constants.c*self.nbins)/(sum(counts)*(max(self.tracked_beam[-1][0])-min(self.tracked_beam[-1][0])))
        print(Ipk)
        return Ipk
    
    # optimizing for minimum energy spread or to specific value at FEL
    def sig_e_opt(self):
        try:
            p_sige = self.pvals.get('p_sige')
            return abs(p_sige - np.sqrt(np.mean(self.tracked_beam[-1][1]**2)))
        except IndexError:
            return np.sqrt(np.mean(self.tracked_beam[-1][1]**2))
        
    def sig_t_e_opt(self):
        return np.sqrt(self.sig_e_opt()**2 + self.sig_t_opt()**2)