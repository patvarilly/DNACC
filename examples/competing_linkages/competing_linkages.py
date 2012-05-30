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
# Python script to produce all the figures in Bortolo's competing
# interactions paper:
#
# B.M. Mognetti, M.E. Leunissen and D. Frenkel, Soft Matter 8, 2213 (2012),
#   doi: 10.1039/c2sm06635a

import numpy as np
from math import sqrt
import subprocess
import operator

import dnacc
from dnacc.units import nm

# Set up basic system
plates = dnacc.PlatesMeanField()
L = 20 * nm
plates.set_tether_type_prototype(L=L, sigma=0.0)
ALPHA = plates.add_tether_type(plate='lower', sticky_end='alpha')
BETA = plates.add_tether_type(plate='lower', sticky_end='beta')
ALPHA_P = plates.add_tether_type(plate='upper', sticky_end='alphap')
BETA_P = plates.add_tether_type(plate='upper', sticky_end='betap')


# A few useful utility methods
def reset_plates(plates):
    for t in plates.tether_types:
        t['sigma'] = 0.0
        t['L'] = L
    plates.beta_DeltaG0.clear()


def set_competing_interactions(plates, beta_DeltaG0a, beta_DeltaG0b):
    plates.beta_DeltaG0['alpha', 'alphap'] = beta_DeltaG0a
    plates.beta_DeltaG0['alpha', 'betap'] = beta_DeltaG0b
    plates.beta_DeltaG0['beta', 'alphap'] = beta_DeltaG0b

# Sample interactions potentials at every 0.05 * L
hArr = np.linspace(0.05 * L, 2.00 * L, 40)


# Figure 2
# ========
def figure2():
    reset_plates(plates)
    S = 0.75 * sqrt(2.0)
    sigma = 1 / (S * L) ** 2
    for t in plates.tether_types:
        t['sigma'] = sigma

    plates.separation = L

    with open('fig2.txt', 'w') as f:
        f.write('# betaDeltaG0a (kT)\t' 'n_alpha / N\t' 'n_beta / N\n')

        for beta_DeltaDeltaG in (3, 5, 8, 11, 14):
            f.write("# beta_DeltaDeltaG = %g kT\n" % beta_DeltaDeltaG)

            for beta_DeltaG0a in xrange(0, -51, -1):
                set_competing_interactions(plates, beta_DeltaG0a,
                                           beta_DeltaG0a + beta_DeltaDeltaG)
                plates.update()

                f.write('%g\t%g\t%g\n' %
                        (beta_DeltaG0a,
                         plates.sigma_bound[ALPHA, ALPHA_P] / (2 * sigma),
                         (plates.sigma_bound[ALPHA, BETA_P] +
                          plates.sigma_bound[BETA, ALPHA_P]) / (2 * sigma)))

            f.write('\n\n')

    subprocess.call(['gnuplot', 'plot_fig2.gp'])


# Figure 3a
# =========
def figure3a():
    single = dnacc.PlatesMeanField()
    single.set_tether_type_prototype(L=L, sigma=0)
    ALPHA = single.add_tether_type(plate='lower', sticky_end='alpha')
    ALPHA_P = single.add_tether_type(plate='upper', sticky_end='alphap')

    with open('fig3a.txt', 'w') as f:
        f.write('# betaDeltaG0a (kT)\t' 'n_alpha / N (MF)\t'
                'n_alpha / N (SCMF)\n')

        for S in (1.06, 0.75, 0.53):
            f.write('# S = %.1f L_alpha\n' % S)
            sigma = 1 / (S * L) ** 2
            for t in single.tether_types:
                t['sigma'] = sigma

            for h in (1, 1.5):
                f.write('# h = %.1f L_alpha\n' % h)
                single.separation = h * L

                for beta_DeltaG0a in xrange(0, -31, -1):
                    single.beta_DeltaG0['alpha', 'alphap'] = beta_DeltaG0a
                    single.update()

                    f.write('%g\t%g\t%g\n' %
                            (beta_DeltaG0a,
                             # Mean Field (MF)  [not really]
                             single.sigma_bound[ALPHA, ALPHA_P] / sigma,
                             # Self-consistent mean field
                             single.sigma_bound[ALPHA, ALPHA_P] / sigma))

                f.write('\n\n')

    subprocess.call(['gnuplot', 'plot_fig3a.gp'])


# Figure 3b
# =========
def figure3b():
    reset_plates(plates)
    S = 0.75
    sigma = 1 / (S * L) ** 2
    for t in plates.tether_types:
        t['sigma'] = sigma

    plates.separation = L

    with open('fig3b.txt', 'w') as f:
        f.write('# betaDeltaG0a (kT)\t' 'n_alpha / N\t' 'n_beta / N\n')

        beta_DeltaDeltaG = 8
        for beta_DeltaG0a in xrange(0, -41, -1):
            set_competing_interactions(plates, beta_DeltaG0a,
                                       beta_DeltaG0a + beta_DeltaDeltaG)
            plates.update()

            f.write('%g\t%g\t%g\n' %
                    (beta_DeltaG0a,
                     plates.sigma_bound[ALPHA, ALPHA_P] / (2 * sigma),
                     (plates.sigma_bound[ALPHA, BETA_P] +
                      plates.sigma_bound[BETA, ALPHA_P]) / (2 * sigma)))

    subprocess.call(['gnuplot', 'plot_fig3b.gp'])


# Figure 4a
# =========
def figure4a():
    reset_plates(plates)
    S = 0.75
    ts = plates.tether_types
    ts[ALPHA]['sigma'] = ts[ALPHA_P]['sigma'] = 0.3 / (S * L) ** 2
    ts[BETA]['sigma'] = ts[BETA_P]['sigma'] = 0.7 / (S * L) ** 2

    with open('fig4a.txt', 'w') as f:
        f.write('# betaDeltaG0b (kT)\t' 'F_min (kT/L^2)\n')

        for beta_DeltaDeltaG in (-1000, 8, 5, 3):
            f.write("# beta_DeltaDeltaG = %g kT\n" % beta_DeltaDeltaG)

            for beta_DeltaG0b in xrange(-20, 9):
                set_competing_interactions(plates,
                                           beta_DeltaG0b - beta_DeltaDeltaG,
                                           beta_DeltaG0b)

                f.write('%g\t%g\n' %
                        (beta_DeltaG0b,
                         min((plates.at(h).free_energy_density
                              for h in hArr)) / (1 / L ** 2)))

            f.write('\n\n')

    subprocess.call(['gnuplot', 'plot_fig4a.gp'])


# Figure 4b
# =========
def figure4b():
    reset_plates(plates)
    S = 0.75
    ts = plates.tether_types
    ts[ALPHA]['sigma'] = ts[ALPHA_P]['sigma'] = 0.3 / (S * L) ** 2
    ts[BETA]['sigma'] = ts[BETA_P]['sigma'] = 0.7 / (S * L) ** 2

    with open('fig4b.txt', 'w') as f:
        f.write('# f (Fraction of hybridised linkages)\t'
                'F_min (kT/L^2)\n')

        for beta_DeltaDeltaG in (+1000, -1000, 8, 5, 3):
            f.write("# beta_DeltaDeltaG = %g kT\n" % beta_DeltaDeltaG)

            for beta_DeltaG0b in xrange(-24, 5):
                beta_DeltaG0a = beta_DeltaG0b - beta_DeltaDeltaG
                if beta_DeltaDeltaG == +1000:
                    set_competing_interactions(plates,
                                               beta_DeltaG0b,
                                               +1000)
                else:
                    set_competing_interactions(plates,
                                               beta_DeltaG0a,
                                               beta_DeltaG0b)

                hAtMin, minF = min(((h, plates.at(h).free_energy_density)
                                    for h in hArr),
                                   key=operator.itemgetter(1))

                plates.at(hAtMin)

                if beta_DeltaDeltaG == +1000:
                    maxFract = plates.tether_types[ALPHA]['sigma']
                else:
                    maxFract = 2 * min(plates.tether_types[ALPHA]['sigma'],
                                       plates.tether_types[BETA]['sigma'])

                fract = sum(plates.sigma_bound[x]
                            for x in [(ALPHA, ALPHA_P), (ALPHA, BETA_P),
                                      (BETA, ALPHA_P)]) / maxFract

                f.write('%g\t%g\n' % (fract, minF / (1 / L ** 2)))

            f.write('\n\n')

    subprocess.call(['gnuplot', 'plot_fig4b.gp'])


# Figure 5a
# =========
#
# Here, the mean field theory seems to work poorly due to the importance
# of coverage fluctuations

def figure5a():
    reset_plates(plates)
    S = 0.75
    sigma = 0.5 / (S * L) ** 2
    ts = plates.tether_types
    ts[ALPHA]['L'] = ts[ALPHA_P]['L'] = 0.3 * L
    ts[BETA]['L'] = ts[BETA_P]['L'] = 1.7 * L
    for t in plates.tether_types:
        t['sigma'] = sigma

    plates.separation = 0.3 * L

    with open('fig5a.txt', 'w') as f:
        f.write('# betaDeltaG0a (kT)\t' 'n_alpha / N\t' 'n_beta / N\n')

        for beta_DeltaDeltaG in (3, 5, 8, 11, 14):
            f.write("# beta_DeltaDeltaG = %g\n" % beta_DeltaDeltaG)

            for beta_DeltaG0a in xrange(-40, 1):
                set_competing_interactions(plates, beta_DeltaG0a,
                                           beta_DeltaG0a + beta_DeltaDeltaG)
                plates.update()

                f.write('%g\t%g\t%g\n' %
                        (beta_DeltaG0a,
                         plates.sigma_bound[ALPHA, ALPHA_P] / (2 * sigma),
                         (plates.sigma_bound[ALPHA, BETA_P] +
                          plates.sigma_bound[BETA, ALPHA_P]) / (2 * sigma)))

            f.write('\n\n')

    subprocess.call(['gnuplot', 'plot_fig5a.gp'])


# Figure 5b
# =========
def figure5b():
    reset_plates(plates)
    S = 0.75
    sigma = 0.5 / (S * L) ** 2
    ts = plates.tether_types
    ts[ALPHA]['L'] = ts[ALPHA_P]['L'] = 0.3 * L
    ts[BETA]['L'] = ts[BETA_P]['L'] = 1.7 * L
    for t in plates.tether_types:
        t['sigma'] = sigma

    with open('fig5b.txt', 'w') as f:
        f.write('# h/L (L = 20 nm)\t' 'F (kT / L^2)\n')

        beta_DeltaDeltaG = 8
        for beta_DeltaG0a in (-5.8, -8.7, -11.6, -14.5, -17.4, -20.3):
            f.write("# beta_DeltaG0a = %g\n" % beta_DeltaG0a)

            set_competing_interactions(plates, beta_DeltaG0a,
                                       beta_DeltaG0a + beta_DeltaDeltaG)

            VArr = [plates.at(h).free_energy_density for h in hArr]
            for (h, V) in zip(hArr, VArr):
                f.write('%g\t%g\n' % (h / L, V / (1 / L ** 2)))

            f.write('\n\n')

    subprocess.call(['gnuplot', 'plot_fig5b.gp'])

# Main module
figure2()
figure3a()
figure3b()
figure4a()
figure4b()
figure5a()
figure5b()
