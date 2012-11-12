# Copyright 2012 Patrick Varilly, Stefano Angioletti-Uberti
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#!/usr/bin/env python
# Python script to produce plate-plate and sphere-sphere potentials
#  for DNACCs with single-stranded DNA tethers.
#
# The tethers are modeled as freely jointed chains of 8 segments, each
#  of length 5 nm.  Bortolo Mognetti has calculated the integrated
#  Boltzmann factors owing to configurational entropy, and stored these
#  in the accompanying files interN.dat, intraRed.dat and 4P.dat
#
# These calculations were done to prepare a response to the following paper:
#
# W.B. Rogers and J.C. Crocker, Proc. Natl. Acad. Sci. USA 108, 15687 (2011),
#   doi: 10.1073/pnas.1109853108

import numpy as np
from math import pi
import subprocess
import scipy.interpolate

import dnacc
from dnacc.units import nm


# Set up TetherStatistics object for pre-computed ssDNA Boltzmann factors
class ssDNAStatistics(object):

    # Read in integrated Boltzmann factors
    l_Kuhn = 5 * nm

    raw = np.loadtxt('interN.dat')
    interp_bridge = scipy.interpolate.interp1d(
        raw[:, 0] * l_Kuhn, raw[:, 1] * l_Kuhn ** 2,
        bounds_error=False, fill_value=0.0)

    raw = np.loadtxt('intraRed.dat')
    interp_loop = scipy.interpolate.interp1d(
        raw[:, 0] * l_Kuhn, raw[:, 1] * l_Kuhn ** 2,
        bounds_error=False, fill_value=raw[-1, 1] * l_Kuhn ** 2)

    raw = np.loadtxt('4P.dat')
    interp_exclude = scipy.interpolate.interp1d(
        raw[:, 0] * 1 * nm, np.exp(-raw[:, 1]),
        bounds_error=False, fill_value=1.0)

    # These methods get called by dnacc code
    @classmethod
    def calc_boltz_binding_cnf_bridge(cls, system,
                                      type_info_i, type_info_j):
        return (float(cls.interp_bridge(system.separation)) *
                float(cls.interp_exclude(system.separation)) ** 2)

    @classmethod
    def calc_boltz_binding_cnf_loop(cls, system,
                                    type_info_i, type_info_j):
        return (float(cls.interp_loop(system.separation)) *
                float(cls.interp_exclude(system.separation)) ** 2)

    @classmethod
    def calc_boltz_exclusion(cls, system, type_info_i):
        return float(cls.interp_exclude(system.separation))

    @classmethod
    def check_system(cls, system):
        if system.separation <= 0:
            raise ValueError("Invalid plate separation")


# Set up basic system
plates = dnacc.PlatesMeanField(ssDNAStatistics)
R = 550 * nm
area = 4 * pi * R ** 2
ALPHA = plates.add_tether_type(
    plate='lower', sigma=4800.0 / area, sticky_end='alpha')
ALPHA_P = plates.add_tether_type(
    plate='upper', sigma=4200.0 / area, sticky_end='alphap')

# Solution hybridisation free energy in kT given by -(c1 / (zr + TinC) - c2)
c1 = 24070.0
c2 = 70.2964
zr = 273.15

# The distances that we'll sample
hArr = np.arange(5 * nm, 81 * nm, 1 * nm)

# Look at various temperatures
for T in (30.5, 32.0, 33.0, 35.0, 36.0):
    beta_DeltaG0 = -(c1 / (zr + T) - c2)
    plates.beta_DeltaG0['alpha', 'alphap'] = beta_DeltaG0

    # Calculate plate potential
    betaFPlate = [plates.at(h).free_energy_density for h in hArr]
    betaFRepPlate = [plates.at(h).rep_free_energy_density for h in hArr]

    # Print out plate potential
    with open('plates-A_B-T%.1f-G%.1f.txt' %
              (T, beta_DeltaG0), 'w') as f:
        f.write('# h (nm)\t' 'F_rep (kT/nm^2)\t' 'F_att (kT/nm^2)\t'
                'F_plate (kT/nm^2)\n')

        for h, betaF, betaFRep in zip(hArr, betaFPlate, betaFRepPlate):
            betaFAtt = betaF - betaFRep
            f.write('%g\t%g\t%g\t%g\n' %
                    (h / nm,
                     betaFRep / (1 / nm ** 2),
                     betaFAtt / (1 / nm ** 2),
                     betaF / (1 / nm ** 2)))

    # Now compute a similar plate potential using the Poisson approximation
    badBetaFPlate = [plates.rep_free_energy_density -
                     plates.at(h).sigma_bound[ALPHA, ALPHA_P] for h in hArr]

    with open('bad-plates-A_B-T%.1f-G%.1f.txt' %
              (T, beta_DeltaG0), 'w') as f:
        f.write('# h (nm)\t' 'F_rep (kT/nm^2)\t' 'F_att (kT/nm^2)\t'
                'F_plate (kT/nm^2)\n')

        for h, betaF, betaFRep in zip(hArr, badBetaFPlate, betaFRepPlate):
            betaFAtt = betaF - betaFRep
            f.write('%g\t%g\t%g\t%g\n' %
                    (h / nm,
                     betaFRep / (1 / nm ** 2),
                     betaFAtt / (1 / nm ** 2),
                     betaF / (1 / nm ** 2)))

    # Now for sphere-sphere potentials
    betaFSpheres = dnacc.calc_spheres_potential(hArr, betaFPlate, R)
    betaFRepSpheres = dnacc.calc_spheres_potential(hArr, betaFRepPlate, R)

    with open('spheres-A_B-T%.1f-G%.1f.txt' %
              (T, beta_DeltaG0), 'w') as f:
        f.write('# h (nm)\t' 'F_rep (kT)\t' 'F_att (kT)\t'
                'F_spheres (kT)\n')
        for h, betaF, betaFRep in zip(hArr, betaFSpheres, betaFRepSpheres):
            betaFAtt = betaF - betaFRep
            f.write('%g\t%g\t%g\t%g\n' %
                    (h / nm, betaFRep, betaFAtt, betaF))

    # Same, but with the Poisson approximation
    badBetaFSpheres = dnacc.calc_spheres_potential(hArr, badBetaFPlate, R)

    with open('bad-spheres-A_B-T%.1f-G%.1f.txt' %
              (T, beta_DeltaG0), 'w') as f:
        f.write('# h (nm)\t' 'F_rep (kT)\t' 'F_att (kT)\t'
                'F_spheres (kT)\n')
        for h, betaF, betaFRep in zip(hArr, badBetaFSpheres, betaFRepSpheres):
            betaFAtt = betaF - betaFRep
            f.write('%g\t%g\t%g\t%g\n' %
                    (h / nm, betaFRep, betaFAtt, betaF))

    # Now convolve with a Gaussian
    # ============================

    # Set up Gaussian kernel and blurring function
    Ncn = 201
    blurSigma = 3 * nm
    filtX = np.linspace(-3 * blurSigma, +3 * blurSigma, Ncn)
    filtDx = filtX[1] - filtX[0]
    filtX += 0.5 * filtDx

    filtY = np.exp(-filtX ** 2 / (2 * blurSigma ** 2))
    filtY /= np.sum(filtY)

    resampHArr = np.arange(hArr[0], hArr[-1], filtDx)

    def blur(origBetaF):
        interpBetaF = scipy.interpolate.interp1d(hArr, origBetaF)
        resampBetaF = interpBetaF(resampHArr)
        resampExpMinusBetaF = np.exp(-resampBetaF)

        paddedInput = np.concatenate((np.zeros(Ncn / 2), resampExpMinusBetaF,
                                      np.ones(Ncn / 2)))
        blurredExpMinusBetaF = np.convolve(paddedInput, filtY, mode='valid')
        blurredBetaF = -np.log(blurredExpMinusBetaF)

        return blurredBetaF

    # Blur good sphere-sphere potentials
    blurredBetaFSpheres = blur(betaFSpheres)
    blurredBetaFRepSpheres = blur(betaFRepSpheres)
    with open('blurred-spheres-A_B-T%.1f-G%.1f.txt' %
              (T, beta_DeltaG0), 'w') as f:
        f.write('# h (nm)\t' 'blurred F_rep (kT)\t' 'blurred F_att (kT)\t'
                'blurred F_spheres (kT)\n')
        for h, betaF, betaFRep in zip(resampHArr, blurredBetaFSpheres,
                                      blurredBetaFRepSpheres):
            betaFAtt = betaF - betaFRep
            f.write('%g\t%g\t%g\t%g\n' %
                    (h / nm, betaFRep, betaFAtt, betaF))

    # Blur "bad" (i.e., Poisson approx) sphere-sphere potentials
    blurredBadBetaFSpheres = blur(badBetaFSpheres)
    with open('blurred-bad-spheres-A_B-T%.1f-G%.1f.txt' %
              (T, beta_DeltaG0), 'w') as f:
        f.write('# h (nm)\t' 'blurred F_rep (kT)\t' 'blurred F_att (kT)\t'
                'blurred F_spheres (kT)\n')
        for h, betaF, betaFRep in zip(resampHArr, blurredBadBetaFSpheres,
                                      blurredBetaFRepSpheres):
            betaFAtt = betaF - betaFRep
            f.write('%g\t%g\t%g\t%g\n' %
                    (h / nm, betaFRep, betaFAtt, betaF))
