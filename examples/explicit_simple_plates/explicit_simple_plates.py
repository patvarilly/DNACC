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
# Python script to compare explicit-tether potentials between two simple
#  plates against mean-field results (and included Monte Carlo simulations)
#
# Note: The Monte Carlo results were calculated using the same grafting
# points by an external code.  The results are in mc_results_h*.dat.  This
# file also contains results from the explicit-tether self-consistent theory,
# but the fixed-point iteration was performed to a looser tolerance (hence
# the slight numerical differences)

import numpy as np
from math import pi, sqrt
import subprocess
import scipy.interpolate

import dnacc
from dnacc.units import nm

# Basic plate properties
L = 20 * nm
S = 1.0
sigma = 1 / (S * L) ** 2

num_tethers = 400
boxL = sqrt(num_tethers / sigma)


def setup_explicit():
    global explicit_plates

    # Generate plate grafting points
    #grafting_pts = boxL * np.random.random_sample( (num_tethers,2) )
    grafting_pts = np.loadtxt('plate.dat', skiprows=4)

    # Set up system
    explicit_plates = dnacc.Plates(boxL, boxL, periodic=True)
    explicit_plates.set_tether_prototype(L=L)

    cur_plate, next_plate = 'lower', 'upper'
    cur_end, next_end = 'alpha', 'alpha_p'

    for pt in grafting_pts:
        explicit_plates.add_tether(plate=cur_plate,
                                   sticky_end=cur_end,
                                   pos=pt)
        cur_plate, next_plate = next_plate, cur_plate
        cur_end, next_end = next_end, cur_end


## First, use explicit-tether system to calculate <n_bonds>
def run_explicit_nbonds():
    for h in np.linspace(0.1 * L, 2.0 * L, 20):

        print 'h = %g L' % (h / L)
        explicit_plates.separation = h
        explicit_plates.beta_DeltaG0['alpha', 'alpha_p'] = 10
        explicit_plates.update()

        with open('results_h%.1f.txt' % (h / L), 'w') as f:

            for beta_DeltaG0 in np.arange(10.0, -20.1, -1):
                print '    beta_DeltaG0 = %g' % beta_DeltaG0
                explicit_plates.beta_DeltaG0['alpha', 'alpha_p'] = \
                    beta_DeltaG0
                explicit_plates.update(DeltaG0_only=True)

                f.write('%g\t%g\n' % (beta_DeltaG0,
                                      explicit_plates.avg_num_bonds))


## Then, use explicit-tether system to calculate interaction potentials
def calc_explicit_potential(beta_DeltaG0):
    print 'Calculating explicit-tether potential at G_0=%g kT' % beta_DeltaG0
    print ' [this can take a while...]'

    explicit_plates.beta_DeltaG0['alpha', 'alpha_p'] = beta_DeltaG0

    with open('pot_DG%.2f.txt' % beta_DeltaG0, 'w') as f:
        f.write('# h/L\t' 'beta F_att^SC\t' 'beta F^SC\n')

        for h in np.linspace(0.1 * L, 2.0 * L, 20):
            print "    h = %g L" % (h / L)
            explicit_plates.at(h)

            betaF = explicit_plates.free_energy
            betaFRep = explicit_plates.rep_free_energy
            betaFAtt = betaF - betaFRep
            f.write('%g\t%g\t%g\n' % (h / L, betaFAtt / boxL ** 2,
                                      betaF / boxL ** 2))


# Now set up equivalent mean-field system
def setup_mean_field():
    global meanfield_plates

    meanfield_plates = dnacc.PlatesMeanField()
    meanfield_plates.set_tether_type_prototype(
        L=L, sigma=num_tethers / boxL ** 2)
    meanfield_plates.add_tether_type(plate='lower', sticky_end='alpha')
    meanfield_plates.add_tether_type(plate='upper', sticky_end='alphap')


# And calculate interaction potentials with it
def calc_meanfield_potential(beta_DeltaG0):
    print 'Calculating mean-field potential at G_0=%g kT' % beta_DeltaG0

    meanfield_plates.beta_DeltaG0['alpha', 'alphap'] = beta_DeltaG0

    with open('mf_pot_DG%.2f.txt' % beta_DeltaG0, 'w') as f:
        f.write('# h/L\t' 'beta F_att^SC\t' 'beta F^SC\n')

        for h in np.linspace(0.1 * L, 2.0 * L, 20):
            meanfield_plates.at(h)

            betaF = meanfield_plates.free_energy_density
            betaFRep = meanfield_plates.rep_free_energy_density
            betaFAtt = betaF - betaFRep
            f.write('%g\t%g\t%g\n' % (h / L, betaFAtt, betaF))


# Main module
# ===========

setup_explicit()
run_explicit_nbonds()
calc_explicit_potential(-5.0)
calc_explicit_potential(-10.0)
calc_explicit_potential(-15.0)

setup_mean_field()
calc_meanfield_potential(-5.0)
calc_meanfield_potential(-10.0)
calc_meanfield_potential(-15.0)
