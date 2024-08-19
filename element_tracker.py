#  elements 
from scipy import constants 
import numpy as np
import munch
from itertools import chain

import optimisation_functions as ofunc

class Element_Track:
    Cgamma = (constants.e**2)/(3*constants.epsilon_0*(constants.m_e*(constants.c**2))**4)

    def __init__(self, beam, Q, beamline_dict):
        # beam in bunch duration [], energy spread [], central energy [MeV] array
        self.init_beam = np.array(beam)
        self.beam = np.array(beam)
        self.Q = Q
        # the beamline dictionary which defines all of the parameters
        self.beamline_dict = beamline_dict
        self.sigt_vals = []

    # accelerating cavity - voltage in megavolts, phase in degrees
    def acceleration(self, frequency, voltage, phase):
        # modify the passed beam
        beam_mod = []
        beam_mod.append(self.beam[0])
        beam_mod.append((((1 + self.beam[1])*self.beam[2] + constants.e*voltage*np.cos(phase - ((2*np.pi*frequency)/constants.c)*self.beam[0]))/(self.beam[2] + constants.e*voltage*np.cos(phase))) - 1)
        beam_mod.append(self.beam[2] + constants.elementary_charge*voltage*np.cos(phase))
        # update initial definition with modified beam
        self.beam = beam_mod
        # return the modified version
        return np.array(beam_mod)

    # compression - like a chicane or arc - give R56, T566, U5666 or calculate if not given
    def compression(self, RC, *args):
        # modify the passed beam
        beam_mod = []
        if len(args) != 0:
            TC = args[0]
            UC = args[1]
            beam_mod.append(self.beam[0] + RC*self.beam[1] + TC*self.beam[1]**2 + UC*self.beam[1]**3)
        else:
            # 4-dipole chicane assumption
            beam_mod.append(self.beam[0] + RC*self.beam[1] + -(3/2)*RC*self.beam[1]**2 + 2*RC*self.beam[1]**3)
        beam_mod.append(self.beam[1])
        beam_mod.append(self.beam[2])
        # update initial beam definition
        self.beam = beam_mod
        # return the modified version
        return beam_mod 
    
    # synchrotron radiation - each value in the array has a different value applied
    #def synchrotron_radiation(self, Rbend):
    #    isrloss = (self.Cgamma/(2*np.pi))*(np.mean(self.beam[2]**4))/Rbend
    #    print(self.Cgamma)
    #    print(isrloss)
    #    print(np.mean(self.beam[2]))
    
    # takes the defined beamline + converts it to elements
    # 
    def beamline_definition(self):
        beamline_definition = munch.Munch({})
        for element in self.beamline_dict:
            if getattr(self.beamline_dict, element).type == 'RF':
                setattr(beamline_definition, element, self.acceleration(getattr(self.beamline_dict, element).f, getattr(self.beamline_dict, element).V, getattr(self.beamline_dict, element).phi))
            elif getattr(self.beamline_dict, element).type == 'BC':
                # handlina as compression can take R56, T566, U566 or just R56
                if len(getattr(self.beamline_dict, element)) == 2:
                    setattr(beamline_definition, element, self.compression(getattr(self.beamline_dict, element).R56))
                else:
                    setattr(beamline_definition, element, self.compression(getattr(self.beamline_dict, element).R56, getattr(self.beamline_dict, element).T566, getattr(self.beamline_dict, element).U5666))
        return beamline_definition
    
    def track_1D(self):
        track = []
        beamline_definition = self.beamline_definition()
        for element in list(self.beamline_dict.keys()):
            beam = getattr(beamline_definition, element)
            track.append(beam)
        return np.array(track)
    
    # arguments are passed as [element][parameter] list
    # runs tracking
    # reproduces the parameter
    def opt_func(self, args, type, *pvals):
        j = 0
        for element in list(self.beamline_dict.keys()):
            for parameter in list(getattr(self.beamline_dict, element).keys())[1:]:
                setattr(getattr(self.beamline_dict, element), parameter, args[j])
                j+=1
        #print(np.sqrt(np.mean(tracked_beam[len(tracked_beam)-1][0]**2))/(constants.c*constants.femto))
        # this next step can be called as a function - have defined optimisatons I, sig_e, sig_t + combined functions (adding in quadrature)
        ofun = ofunc.Optimisation_Functions(self.track_1D(), *pvals)

        if type == "sig_t":
            # reset the beam
            self.beam = self.init_beam
            return ofun.sig_t_opt()
            #return ofun.sig_t_opt()
        elif type == "sig_e":
            self.beam = self.init_beam
            return ofun.sig_e_opt()
        elif type =='Ipk':
            self.beam = self.init_beam
            return ofun.Ipk_opt(self.Q)
        elif type == "sig_t_e" or "sig_e_t":
            self.beam = self.init_beam
            return ofun.sig_t_e_opt()
        else:
            print("please enter a correct optimisation type")

    
    def initial_guess(self):
        params = []
        for element in list(self.beamline_dict.keys()):
            params.append(list(self.beamline_dict[element].values())[1:])
        initial_guess = list(chain.from_iterable(params))
        return initial_guess
