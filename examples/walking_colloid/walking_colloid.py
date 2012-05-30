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
# Python script to predict properties of DNACCs walking on DNA coverage
# gradients
#
# Written by Patrick Varilly, 7 Nov 2011
# Adapted to new dnacc module, Feb 2012

import numpy as np
import subprocess
import operator

import dnacc
from dnacc.units import nm

# Set up basic system
plates = dnacc.PlatesMeanField()
L = 20 * nm
plates.set_tether_type_prototype(L=L, sigma=0.0)

ALPHA = plates.add_tether_type(plate='walker', sticky_end='alpha')
BETA_1 = plates.add_tether_type(plate='surface', sticky_end='beta1')
BETA_2 = plates.add_tether_type(plate='surface', sticky_end='beta2')

# The rolling walker has this radius
R = 500.0 * nm

# Look at various overall binding energies
ts = plates.tether_types


def do_it(beta_DeltaG0Mid):

    print 'Working on beta_DeltaG0Mid = %g' % beta_DeltaG0Mid

    # and various energy gaps
    for beta_Delta in xrange(0, 10):

        plates.beta_DeltaG0['alpha', 'beta1'] = \
            beta_DeltaG0Mid - 0.5 * beta_Delta
        plates.beta_DeltaG0['alpha', 'beta2'] = \
            beta_DeltaG0Mid + 0.5 * beta_Delta

        # and various total strand coverage
        for S in 0.75, 0.25:

            sigma = 1 / (S * L) ** 2
            ts[ALPHA]['sigma'] = sigma * 0.5

            with open('walk-S%0.2f-G0Mid%.1f-delta%.1f.dat' %
                      (S, beta_DeltaG0Mid, beta_Delta), 'w') as f:

                f.write('c\t' 'F_rep (kT)\n')
                offset = 0
                for c in np.linspace(1.0, 0.0, 21):

                    ts[BETA_1]['sigma'] = c * sigma
                    ts[BETA_2]['sigma'] = (1 - c) * sigma

                    # Calculate effective plate-plate potential
                    hArr = np.linspace(0.05 * L, 2.00 * L, 40)
                    VArr = [plates.at(h).free_energy_density
                            for h in hArr]

                    # and the resulting plate-sphere potential
                    Rplt = 1000.0 * R
                    betaF = dnacc.calc_spheres_potential(
                        hArr, VArr, R, Rplt)

                    # Find its minimum
                    minH, minBetaF = min(zip(hArr, betaF),
                                         key=operator.itemgetter(1))

                    if offset == 0:
                        offset = -minBetaF

                    f.write('%.2f\t%.4f\t%.4f\n' %
                            (c, minBetaF + offset, minH / nm))

for beta_DeltaG0Mid in xrange(-20, 1):
    do_it(beta_DeltaG0Mid)
