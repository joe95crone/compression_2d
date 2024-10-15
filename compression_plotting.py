# plotting functions - shoult only take in the initial and tracked beam phase spaces and should plot from that

import numpy as np
from scipy import constants
from scipy import stats
from matplotlib import pyplot as plt

# may be a few more arguments for peak charge calculation etc.
# one tracking run defines one set of plotting fucntions
class Track_Plot:
    def __init__(self, beamline_dict, init_beam, tracked_beam, label, Q, Npart, nbins = 100, save_plots=False):
        self.beamline = ["Initial"] + list(beamline_dict.keys())
        self.init_beam = init_beam
        self.tracked_beam = tracked_beam
        
        self.label = str(label)
        self.Q = Q
        self.Npart = Npart
        self.nbins = nbins

        self.save_plots = save_plots

    def beam_energy_plt(self):
        beam_energy = []
        for ele in range(len(self.tracked_beam)):
            beam_energy.append(np.mean(self.tracked_beam[ele][2])/(constants.elementary_charge*constants.mega))
        beam_energy = np.insert(beam_energy, 0, np.mean(self.init_beam[2])/(constants.elementary_charge*constants.mega), axis=0)
        plt.figure()
        plt.plot(self.beamline, beam_energy, label=self.label)
        plt.title(self.label + " Beamline Energy")
        plt.xlabel("Element")
        plt.ylabel("Beam Energy [MeV]")
        plt.legend()
        if self.save_plots == True:
            filename = 'FIGS\\' + self.label + '_beamline_energy' + '.png'
            plt.savefig(filename)
        plt.draw()
        return beam_energy

    def beam_sigt_plt(self):
        beam_sigt = []
        for ele in range(len(self.tracked_beam)):
            beam_sigt.append(np.sqrt(np.mean(self.tracked_beam[ele][0]**2))/(constants.c*constants.femto))
        beam_sigt = np.insert(beam_sigt, 0, np.sqrt(np.mean(self.init_beam[0]**2))/(constants.c*constants.femto), axis=0)
        plt.figure()
        plt.plot(self.beamline, beam_sigt, label=self.label)
        plt.title(self.label + " Beamline Bunch Duration")
        plt.xlabel("Element")
        plt.ylabel("Bunch Duration [fs]")
        plt.legend()
        if self.save_plots == True:
            filename = 'FIGS\\' + self.label + '_beamline_bunch_duration' + '.png'
            plt.savefig(filename)
        plt.draw()
        return beam_sigt

    def beam_sige_plt(self):
        beam_sige = []
        for ele in range(len(self.tracked_beam)):
            beam_sige.append(np.sqrt(np.mean(self.tracked_beam[ele][1]**2))*100)
        beam_sige = np.insert(beam_sige, 0, np.sqrt(np.mean(self.init_beam[1]**2)*100), axis=0)
        plt.figure()
        plt.plot(self.beamline, beam_sige, label=self.label)
        plt.title(self.label + " Beamline Energy Spread")
        plt.xlabel("Element")
        plt.ylabel("Energy Spread [%]")
        plt.legend()
        if self.save_plots == True:
            filename = 'FIGS\\' + self.label + '_beamline_energy_spread' + '.png'
            plt.savefig(filename)
        plt.draw()
        return beam_sige

    def beam_Ipk_plt(self):
        beam_Ipk = []
        for ele in range(len(self.tracked_beam)):
            counts, bins = np.histogram(self.tracked_beam[ele][0]/(constants.c*constants.femto), self.nbins)
            beam_Ipk.append(max(counts)*self.bin_current(ele))
        init_counts, init_bins = np.histogram(self.init_beam[0]/(constants.c*constants.femto), self.nbins)
        init_Ipk = max(init_counts)*self.bin_current(0)
        beam_Ipk = np.insert(beam_Ipk, 0, init_Ipk, axis=0)
        plt.figure()
        plt.plot(self.beamline, beam_Ipk/constants.kilo, label=self.label)
        plt.title(self.label + " Beamline Peak Current")
        plt.xlabel("Element")
        plt.ylabel("Peak Current [kA]")
        plt.legend()
        if self.save_plots == True:
            filename = 'FIGS\\' + self.label + '_beamline_peak_current' + '.png'
            plt.savefig(filename)
        plt.draw()
        return beam_Ipk

    # calculate the Pearson correlation coefficient throughout the FEL
    def beam_linearity_plt(self):
        pearson_linearity = []
        for ele in range(len(self.tracked_beam)):
            pearson_linearity.append(stats.pearsonr(self.tracked_beam[ele][0], self.tracked_beam[ele][1], alternative = 'two-sided').statistic)
        pearson_linearity = np.insert(pearson_linearity, 0, stats.pearsonr(self.init_beam[0], self.init_beam[1], alternative = 'two-sided').statistic, axis=0)
        plt.figure()
        plt.plot(self.beamline, pearson_linearity, label=self.label)
        plt.title(self.label + " Beamline LPS Linearity")
        plt.xlabel("Element")
        plt.ylabel("Pearson Correlation Coefficient")
        plt.legend()
        if self.save_plots == True:
            filename = 'FIGS\\' + self.label + '_linearity' + '.png'
            plt.savefig(filename)
        plt.draw()
        return pearson_linearity

    def beam_sigt_sige_phase_space(self, position=-1):
        plt.figure()
        if position == 0:
            plt.scatter(self.init_beam[0]/(constants.c*constants.femto), self.init_beam[1], label=self.label)
        elif position == -1:
            plt.scatter(self.tracked_beam[-1][0]/(constants.c*constants.femto), self.tracked_beam[-1][1], label=self.label)
        else:
            plt.scatter(self.tracked_beam[position-1][0]/(constants.c*constants.femto), self.tracked_beam[position-1][1], label=self.label)
        plt.title(self.label + " " + self.beamline[position])
        plt.xlabel("T [fs]")
        plt.ylabel("$\mathregular{\sigma_{E}}$ ")
        plt.legend()
        if self.save_plots == True:
            filename = 'FIGS\\' + self.label + '_sigT_sigE_' + self.beamline[position] + '.png'
            plt.savefig(filename)
        plt.draw()

    def beam_sigz_sige_phase_space(self, position=-1):
        plt.figure()
        if position == 0:
            plt.scatter(self.init_beam[0]/constants.milli, self.init_beam[1], label=self.label)
        elif position == -1:
            plt.scatter(self.tracked_beam[-1][0]/constants.milli, self.tracked_beam[-1][1], label=self.label)
        else:
            plt.scatter(self.tracked_beam[position-1][0]/constants.milli, self.tracked_beam[position-1][1], label=self.label)
        plt.title(self.label + " " + self.beamline[position])
        plt.xlabel("Z [mm]")
        plt.ylabel("$\mathregular{\sigma_{E}}$ ")
        plt.legend()
        if self.save_plots == True:
            filename = 'FIGS\\' + self.label + '_sigT_sigE_' + self.beamline[position] + '.png'
            plt.savefig(filename)
        plt.draw()

    # factor to convert counts to current
    def bin_current(self, position=-1):
        if position == 0:
            return (self.Q*constants.pico*constants.c*self.nbins)/(self.Npart*(max(self.init_beam[0])-min(self.init_beam[0])))
        elif position == -1:
            return (self.Q*constants.pico*constants.c*self.nbins)/(self.Npart*(max(self.tracked_beam[-1][0])-min(self.tracked_beam[-1][0])))
        else:
            return (self.Q*constants.pico*constants.c*self.nbins)/(self.Npart*(max(self.tracked_beam[position-1][0])-min(self.tracked_beam[position-1][0])))
        
    # beam current histogram
    def current_plt(self, position=-1):
        if position == 0:
            counts, bins = np.histogram(self.tracked_beam[0]/(constants.c*constants.femto), self.nbins)
        elif position == -1:
            counts, bins = np.histogram(self.tracked_beam[-1][0]/(constants.c*constants.femto), self.nbins)
        else:
            counts, bins = np.histogram(self.tracked_beam[position-1][0]/(constants.c*constants.femto), self.nbins)
        Ipk = max(counts)*self.bin_current(position)
        plt.figure()
        plt.stairs(counts*self.bin_current(position), bins, label=self.label)
        plt.title("{0} {1} Ipk = {2:.2f} kA".format(self.label, self.beamline[position], Ipk/constants.kilo))
        plt.xlabel("T [fs]")
        plt.ylabel("Current [A]")
        plt.legend()
        if self.save_plots == True:
            filename = 'FIGS\\' + self.label + '_Current_' + self.beamline[position] + '.png'
            plt.savefig(filename)
        plt.draw() 
        return Ipk
        