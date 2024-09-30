# beam generation
import numpy as np
from scipy import constants
import h5py
import os

class HDF5_read:
	def __init__(self):
		# Store an empty list for dataset names
		self.names = []

	def __call__(self, name, h5obj):
		# only h5py datasets have dtype attribute, so we can search on this
		if hasattr(h5obj,'dtype') and not name in self.names:
			self.names += [name]

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
        elif init_type == "file":
            self.LPS_data_tag = ['q','z','pz']
            self.filename = args[0]
            self.path = os.getcwd() + '\\Distributions\\'
            self.file = h5py.File(self.path + self.filename,'r')
            self.HDF5_in = HDF5_read()
            self.file.visititems(self.HDF5_in)

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
    
    # generate a cubic bunch in 1D using a series of dz values
    def cubic_sige_generator(self, dz):
        scale = (self.sig_e/(self.sig_t*constants.c))
        fac = 0.85/(self.sig_e/1e-4)
        return -scale*dz + (1/(fac*scale)**3)*dz**3

    # generate a Gaussian bunch
    def gaussian_bunch(self):
        beam_sigz = np.array(self.gaussian_generator(0, self.sig_t*constants.c, self.cutoff))
        beam_sige = np.array(self.gaussian_generator(0, self.sig_e, self.cutoff))
        beam_E = np.full(self.npart, self.E_mean)
        beam = [beam_sigz, beam_sige, beam_E]
        return beam

    # generate a more realistic bunch with a cubic term in the energy spread
    def cubic_sige_bunch(self):
        beam_sigz = np.array(self.gaussian_generator(0, self.sig_t*constants.c, self.cutoff))
        beam_sige = np.array(self.cubic_sige_generator(beam_sigz))
        beam_E = np.full(self.npart, self.E_mean)
        beam = [beam_sigz, beam_sige, beam_E]
        return beam
    
    # get the bunch data in the standard three column format
    # also sets spome of the parameters of the class
    #! Unit conversion is based on Can's standard OPAL input 
    def get_file_LPS_data(self):
        dat_names = [] 
        for tag in range(len(self.LPS_data_tag)):
            for name in range(len(self.HDF5_in.names)):
                if self.LPS_data_tag[tag] == self.HDF5_in.names[name].split('/')[1]:
                    dat_names.append(self.HDF5_in.names[name])
        HDF5_LPS = np.array([np.array(self.file[dat_names[i]]) for i in range(len(dat_names))])
        # no macroparticles
        self.npart = len(HDF5_LPS[0])
        # charge (pC)
        self.bunch_charge = -sum(HDF5_LPS[0])/constants.pico 
        # SI Emean
        self.Emean = np.mean(HDF5_LPS[2])*constants.value('electron mass energy equivalent in MeV')*constants.elementary_charge*constants.mega
        delta_e = np.array([(HDF5_LPS[2][i]-np.mean(HDF5_LPS[2]))/np.mean(HDF5_LPS[2]) for i in range(len(HDF5_LPS[2]))])
        Emean_beam = np.full(len(HDF5_LPS[2]), self.Emean)
        beam = [np.array(HDF5_LPS[1]), delta_e, Emean_beam]
        return beam

#if __name__ == "__main__":
    # example generation run
    #beam_gen_param = Beam_Generation("param", 250, 1e4, 3, 1, 1e-4, 100)
    #print(beam_gen_param.bunch_charge)
    #Gbeam = beam_gen_param.gaussian_bunch()
    #print(Gbeam)

    # example from file run
    #beam_gen = Beam_Generation("file","#30.h5")
    #beam = beam_gen.get_file_LPS_data()
    #print(beam_gen.bunch_charge)
    #print(beam)