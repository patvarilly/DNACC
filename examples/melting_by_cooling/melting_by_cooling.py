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
# Python script to produce all the figures in Stefano's melting-by-cooling
# paper
#
# S. Angioletti-Uberti, B.M. Mognetti and D. Frenkel, "Reentrant melting
# as a design principle for DNA-coated colloids", Nature Materials (11),
# 518--522 (2012)

import numpy as np
from math import sqrt
import subprocess
import operator

import dnacc
from dnacc.units import nm

# Set up basic system
platesXX = dnacc.PlatesMeanField()
L = 20 * nm

platesXX.set_tether_type_prototype(L=L, sigma=0.0, plate='lower')
XX_LOWER_ALPHA = platesXX.add_tether_type(sticky_end='alpha')
XX_LOWER_BETA = platesXX.add_tether_type(sticky_end='beta')

platesXX.set_tether_type_prototype(L=L, sigma=0.0, plate='upper')
XX_UPPER_ALPHA = platesXX.add_tether_type(sticky_end='alpha')
XX_UPPER_BETA = platesXX.add_tether_type(sticky_end='beta')


platesXXP = dnacc.PlatesMeanField()
L = 20 * nm

platesXXP.set_tether_type_prototype(L=L, sigma=0.0, plate='lower')
XXP_LOWER_ALPHA = platesXXP.add_tether_type(sticky_end='alpha')
XXP_LOWER_BETA = platesXXP.add_tether_type(sticky_end='beta')

platesXXP.set_tether_type_prototype(L=L, sigma=0.0, plate='upper')
XXP_UPPER_ALPHA_P = platesXXP.add_tether_type(sticky_end='alphap')
XXP_UPPER_BETA_P = platesXXP.add_tether_type(sticky_end='betap')


# A few useful utility methods
def reset_plates(plates, S):
    for t in plates.tether_types:
        t['sigma'] = 1 / (S * L) ** 2
    plates.beta_DeltaG0.clear()


def set_competing_interactions(plates, beta_DeltaG0a, beta_DeltaG0b):

    plates.beta_DeltaG0['alpha', 'alphap'] = beta_DeltaG0a
    plates.beta_DeltaG0['alpha', 'beta'] = beta_DeltaG0b
    plates.beta_DeltaG0['alphap', 'betap'] = beta_DeltaG0b

# Sample interactions potentials at every 0.05 * L
hArr = np.linspace(0.05 * L, 2.00 * L, 40)


# Figure 2a
# =========
def figure2a():
    reset_plates(platesXXP, 1.0)
    sigma = platesXXP.tether_types[0]['sigma']

    with open('fig2a.txt', 'w') as f:
        f.write('# DeltaG0strong (kT)\t' 'n_strong / N\t'
                'n_weak / N\t' 'F_min/A (kT / L^2)\n')

        beta_Delta = 5

        for beta_DeltaG0strong in xrange(-40, 7):

            set_competing_interactions(platesXXP,
                                       beta_DeltaG0strong,
                                       beta_DeltaG0strong + beta_Delta)
            platesXXP.at(2.1 * L).set_reference_now()

            VArr = [platesXXP.at(h).free_energy_density for h in hArr]
            hAtMin, minF = min(zip(hArr, VArr), key=operator.itemgetter(1))
            platesXXP.at(hAtMin)

            f.write('%g\t%g\t%g\t%g\n' %
                    (beta_DeltaG0strong,
                     platesXXP.sigma_bound[XXP_LOWER_ALPHA,
                                           XXP_UPPER_ALPHA_P] / sigma,
                     (platesXXP.sigma_bound[XXP_LOWER_ALPHA,
                                            XXP_LOWER_BETA] +
                      platesXXP.sigma_bound[XXP_UPPER_ALPHA_P,
                                            XXP_UPPER_BETA_P]) / (2 * sigma),
                     minF / (1 / L ** 2)))

    subprocess.call(['gnuplot', 'plot_fig2a.gp'])


# Figure 2b
# =========
def figure2b():
    reset_plates(platesXX, 1.0)
    sigma = platesXX.tether_types[0]['sigma']

    with open('fig2b.txt', 'w') as f:
        f.write('# DeltaG0strong (kT)\t'
                'n_weak / N\t' 'F_min/A (kT / L^2)\n')

        beta_Delta = 5

        for beta_DeltaG0strong in xrange(-40, 7):

            set_competing_interactions(platesXX,
                                       beta_DeltaG0strong,
                                       beta_DeltaG0strong + beta_Delta)
            platesXX.at(2.1 * L).set_reference_now()

            VArr = [platesXX.at(h).free_energy_density for h in hArr]
            hAtMin, minF = min(zip(hArr, VArr), key=operator.itemgetter(1))
            platesXX.at(hAtMin)

            f.write('%g\t%g\t%g\n' %
                    (beta_DeltaG0strong,
                     (platesXX.sigma_bound[XX_LOWER_ALPHA, XX_LOWER_BETA] +
                      platesXX.sigma_bound[XX_LOWER_ALPHA, XX_UPPER_BETA] +
                      platesXX.sigma_bound[XX_UPPER_ALPHA, XX_LOWER_BETA] +
                      platesXX.sigma_bound[XX_UPPER_ALPHA, XX_UPPER_BETA]) /
                     (2 * sigma),
                     minF / (1 / L ** 2)))

    subprocess.call(['gnuplot', 'plot_fig2b.gp'])


# Figure 2c
# =========
# Plots of the potentials, not actually a figure in the paper
def figure2c():
    reset_plates(platesXXP, 1.0)
    sigma = platesXXP.tether_types[0]['sigma']

    with open('fig2c.txt', 'w') as f:
        f.write('# h/L (L = 20 nm)\t' 'F/A (kT / L^2)\n')

        beta_Delta = 5

        for beta_DeltaG0strong in xrange(-40, 7):

            f.write('# betaDeltaG0strong = %g\n' % beta_DeltaG0strong)

            set_competing_interactions(platesXXP,
                                       beta_DeltaG0strong,
                                       beta_DeltaG0strong + beta_Delta)
            platesXXP.at(2.1 * L).set_reference_now()

            VArr = [platesXXP.at(h).free_energy_density for h in hArr]
            for h, V in zip(hArr, VArr):
                f.write('%g\t%g\n' % (h / L, V / (1 / L ** 2)))

            f.write('\n\n')

    subprocess.call(['gnuplot', 'plot_fig2c.gp'])


# Figure 3a
# =========
def figure3a():
    with open('fig3a.txt', 'w') as f:
        f.write('# DeltaG0strong (kT)\t'
                'F_min_XXp/A (kT / L^2)\t' 'F_min_XX/A (kT / L^2)\n')

        beta_Delta = 5

        for S in (1.5, 1.0, 0.7, 0.5):
            reset_plates(platesXX, S)
            reset_plates(platesXXP, S)

            for beta_DeltaG0strong in xrange(-30, 3):

                set_competing_interactions(platesXX,
                                           beta_DeltaG0strong,
                                           beta_DeltaG0strong + beta_Delta)
                platesXX.at(2.1 * L).set_reference_now()

                set_competing_interactions(platesXXP,
                                           beta_DeltaG0strong,
                                           beta_DeltaG0strong + beta_Delta)
                platesXXP.at(2.1 * L).set_reference_now()

                minFXX = min(platesXX.at(h).free_energy_density
                             for h in hArr)
                minFXXP = min(platesXXP.at(h).free_energy_density
                              for h in hArr)

                f.write('%g\t%g\t%g\n' %
                        (beta_DeltaG0strong,
                         minFXXP / (1 / L ** 2), minFXX / (1 / L ** 2)))

            f.write('\n\n')

    subprocess.call(['gnuplot', 'plot_fig3a.gp'])


# Figure 3b
# =========
# For this figure, the colloids are a bit more complicated, there
# is an additional type of strand there just for steric repulsion

def figure3b():

    # Add extra tether types (removed at the end)
    XX_LOWER_GAMMA = platesXX.add_tether_type(plate='lower',
                                              sticky_end='GAMMA')
    XX_UPPER_GAMMA = platesXX.add_tether_type(plate='upper',
                                              sticky_end='GAMMA')

    XXP_LOWER_GAMMA = platesXXP.add_tether_type(plate='lower',
                                                sticky_end='GAMMA')
    XXP_UPPER_GAMMA = platesXXP.add_tether_type(plate='upper',
                                                sticky_end='GAMMA')

    # Go
    S = 1.0
    reset_plates(platesXX, S)
    reset_plates(platesXXP, S)

    with open('fig3b.txt', 'w') as f:
        f.write('# DeltaG0strong (kT)\t'
                'F_min_XXp/A (kT / L^2)\t' 'F_min_XX/A (kT / L^2)\n')

        beta_Delta = 5

        for Lincrt in (x * L for x in (1.0, 1.2, 1.5, 1.8)):

            platesXX.tether_types[XX_LOWER_GAMMA]['L'] = Lincrt
            platesXX.tether_types[XX_UPPER_GAMMA]['L'] = Lincrt
            platesXXP.tether_types[XXP_LOWER_GAMMA]['L'] = Lincrt
            platesXXP.tether_types[XXP_UPPER_GAMMA]['L'] = Lincrt

            for beta_DeltaG0strong in xrange(-30, 3):

                set_competing_interactions(platesXX,
                                           beta_DeltaG0strong,
                                           beta_DeltaG0strong + beta_Delta)
                platesXX.at(2.1 * L).set_reference_now()

                set_competing_interactions(platesXXP,
                                           beta_DeltaG0strong,
                                           beta_DeltaG0strong + beta_Delta)
                platesXXP.at(2.1 * L).set_reference_now()

                minFXX = min(platesXX.at(h).free_energy_density
                             for h in hArr)
                minFXXP = min(platesXXP.at(h).free_energy_density
                              for h in hArr)

                f.write('%g\t%g\t%g\n' %
                        (beta_DeltaG0strong,
                         minFXXP / (1 / L ** 2), minFXX / (1 / L ** 2)))

            f.write('\n\n')

    # Remove extra GAMMA tethers
    del platesXX.tether_types[XX_LOWER_GAMMA:XX_UPPER_GAMMA + 1]
    del platesXXP.tether_types[XXP_LOWER_GAMMA:XXP_UPPER_GAMMA + 1]
    reset_plates(platesXX, S)
    reset_plates(platesXXP, S)

    subprocess.call(['gnuplot', 'plot_fig3b.gp'])


# Figure 3c
# =========
def figure3c():
    with open('fig3c.txt', 'w') as f:
        f.write('# DeltaG0strong (kT)\t'
                'F_min_XXp/A (kT / L^2)\t' 'F_min_XX/A (kT / L^2)\n')

        S = 0.5
        for beta_Delta in (1, 3, 5, 7):
            f.write('# beta_Delta = %g\n' % beta_Delta)

            reset_plates(platesXX, S)
            reset_plates(platesXXP, S)

            for beta_DeltaG0strong in xrange(-30, 3):

                set_competing_interactions(platesXX,
                                           beta_DeltaG0strong,
                                           beta_DeltaG0strong + beta_Delta)
                platesXX.at(2.1 * L).set_reference_now()

                set_competing_interactions(platesXXP,
                                           beta_DeltaG0strong,
                                           beta_DeltaG0strong + beta_Delta)
                platesXXP.at(2.1 * L).set_reference_now()

                minFXX = min(platesXX.at(h).free_energy_density
                             for h in hArr)
                minFXXP = min(platesXXP.at(h).free_energy_density
                              for h in hArr)

                f.write('%g\t%g\t%g\n' %
                        (beta_DeltaG0strong,
                         minFXXP / (1 / L ** 2), minFXX / (1 / L ** 2)))

            f.write('\n\n')

    subprocess.call(['gnuplot', 'plot_fig3c.gp'])


# Main module
figure2a()
figure2b()
figure2c()
figure3a()
figure3b()
figure3c()
