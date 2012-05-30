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
# Python script to produce effective pair potentials between plates and
#  spheres, a la:
#
# M.E. Leunissen and D. Frenkel, J. Chem. Phys. 134, 084702 (2011),
#   doi: 10.1063/1.3557794
#
# Note that Mirjam's original expression for the configurational
#  contribution to the binding free energy is off by a factor of
#  sin( gamma_i + gamma_j ).  For details, see
#
# B.M. Mognetti, M.E. Leunissen and D. Frenkel, Soft Matter 8, 2213 (2012),
#   doi: 10.1039/c2sm06635a
#
# Some of the remaining discrepancy is due to using mean-field theory and
# to Mirjam's interpolation formula breaking down when the strongest
# plate-plate interaction energy density exceeds 2 kT / L^2


import dnacc
from dnacc.units import nm
import numpy as np

# Set up a simple mean-field plates system
plates = dnacc.PlatesMeanField()

# Add two tether types, one on each plate, with equal lengths
L = 20 * nm
plates.set_tether_type_prototype(sigma=0, L=L)
ALPHA = plates.add_tether_type(plate='upper', sticky_end='alpha')
ALPHA_P = plates.add_tether_type(plate='lower', sticky_end='alphap')

# Where will we calculate these potentials
hArr = np.linspace(1 * nm, 40 * nm, 40)

# Do various grafting densities
for S in 0.25, 0.75:
    sigma = 1 / (S * L) ** 2
    plates.tether_types[ALPHA]['sigma'] = sigma
    plates.tether_types[ALPHA_P]['sigma'] = sigma

    for betaDeltaG0 in xrange(-7, -2):
        plates.beta_DeltaG0['alpha', 'alphap'] = betaDeltaG0

        betaFPlate = [plates.at(h).free_energy_density for h in hArr]

        # Make plate potential first
        with open('plates-S%0.2f-G%.1f.dat' % (S, betaDeltaG0), 'w') as f:

            f.write('\t'.join(['h / L',
                               "F_rep (kT/L^2)",
                               "F_att (kT/L^2)",
                               "F_plate (kT/L^2)"]) + '\n')

            for h, V in zip(hArr, betaFPlate):
                betaFRep = plates.at(h).rep_free_energy_density
                betaFAtt = V - betaFRep

                f.write('%.7g\t%.7g\t%.7g\t%.7g\n'
                        % (h / L, betaFRep / (1 / L ** 2),
                           betaFAtt / (1 / L ** 2),
                           (betaFRep + betaFAtt) / (1 / L ** 2)))

        # Now sphere potentials
        # (don't decompose into attractive and repulsive contributions,
        # although it's not hard to do)
        for R in 6.7, 25.0:
            betaFSphere = dnacc.calc_spheres_potential(
                hArr, betaFPlate, R * L)
            with open('spheres-R%.1f-S%0.2f-G%.1f.dat'
                      % (R, S, betaDeltaG0), 'w') as f:

                f.write('\t'.join(['h / L',
                                   "[ignore: F_rep (kT)]",
                                   "[ignore: F_att (kT)]",
                                   "F_sphere (kT)"]) + '\n')
                for h, V in zip(hArr, betaFSphere):
                    f.write('%.7g\t%.7g\t%.7g\t%.7g\n'
                            % (h / L, 0, 0, V))
