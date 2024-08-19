# beam generation
import numpy as np
from scipy import constants

class Beam_Generation:
    # want to pass either a filename or an image or a set of parameters
    # init_type specifies the beam creation method
    def __init__(self, init_type, *args):    
        # passing a set of parameters
        if init_type == "param" and len(args) == 6:
            self.bunch_charge = args[0] # pC bunch charge
            self.npart = int(args[1]) # no. particles to generate
            self.cutoff = args[2] # Gaussian cut-off truncations (no. sigmas)
            self.sig_t = args[3]*constants.pico # s rms bunch duration 
            self.sig_e = args[4] # rms relative energy spread 
            self.E_mean = args[5]*constants.elementary_charge # SI mean energy
        
    # Gaussian sampling with a cut off
    def gaussian_generator(self, mean, std, trunc):
        samples = []
        for particle in range(self.npart):
            if std == 0:
                samples.append(mean)
            else:
                lims = [mean - (trunc*std), mean + (trunc*std)]
                sample = np.random.normal(mean, std)
                while lims[0] >= sample >= lims[1]:
                    sample = np.random.normal(mean, std) 
                samples.append(sample)
        return samples
    
    def gaussian_bunch(self):
        beam_sigz = np.array(self.gaussian_generator(0, self.sig_t*constants.c, self.cutoff))
        beam_sige = np.array(self.gaussian_generator(0, self.sig_e, self.cutoff))
        beam_E = np.full(self.npart, self.E_mean)
        beam = [beam_sigz, beam_sige, beam_E]
        return beam

#if __name__ == "__main__":
#    beam_gen = Beam_Generation("param", 250, 1e4, 3, 1, 1e-4, 100)
#    print(beam_gen.bunch_charge)
#    Gbeam = beam_gen.gaussian_bunch()