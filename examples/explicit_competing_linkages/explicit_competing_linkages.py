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
# Python script to compare explicit-tether bond probabilities between two
#  plates with competing linkages against Monte Carlo results
#
# See examples/competingLinkages for the mean-field analogues
#
# Note: The Monte Carlo results were calculated using the same grafting
# points by Bortolo's Fortran program, given in CompetingData/4p.  The
# results are in CompetingData/nX_G, where X = a/b for strong/weak bonds
# (alpha-alpha' vs [alpha-beta' and beta-alpha']) and G is the difference in
# binding strength between strong and weak bonds.  Each file has two
# columns: the value of beta_DeltaG0 for strong bonds and the corresponding
# number of X-type bonds formed relative to the maximum possible number of
# weak bonds.

import numpy as np
from math import pi, sqrt
import subprocess
import scipy.interpolate

import dnacc
from dnacc.units import nm

# Basic plate properties
L = 20 * nm
S = 0.75 * sqrt(2.0)
sigma = 1 / (S * L) ** 2
NAlphas = 500
num_tethers = 4 * NAlphas
boxL = sqrt(NAlphas / sigma)


def setup():
    global plates, ALPHAS, ALPHA_PS, BETAS, BETA_PS

    # Assign grafting points

    #grafting_pts_ALPHA = boxL * np.random.random_sample( (NAlphas,2) )
    #grafting_pts_ALPHA_P = boxL * np.random.random_sample( (NAlphas,2) )
    #grafting_pts_BETA = boxL * np.random.random_sample( (NAlphas,2) )
    #grafting_pts_BETA_P = boxL * np.random.random_sample( (NAlphas,2) )

    all_grafting_pts = L * np.loadtxt('CompetingData/4p', skiprows=1)
    grafting_pts_ALPHA = all_grafting_pts[:NAlphas, 0:2]
    grafting_pts_ALPHA_P = all_grafting_pts[:NAlphas, 2:4]
    grafting_pts_BETA = all_grafting_pts[NAlphas:2 * NAlphas, 0:2]
    grafting_pts_BETA_P = all_grafting_pts[NAlphas:2 * NAlphas, 2:4]

    # Set up system
    plates = dnacc.Plates(boxL, boxL, periodic=True)

    ALPHAS, ALPHA_PS, BETAS, BETA_PS = set(), set(), set(), set()

    plates.set_tether_prototype(plate='lower', L=L, sigma=sigma)
    for pt in grafting_pts_ALPHA:
        ALPHAS.add(plates.add_tether(sticky_end='alpha', pos=pt))
    for pt in grafting_pts_BETA:
        BETAS.add(plates.add_tether(sticky_end='beta', pos=pt))

    plates.set_tether_prototype(plate='upper', L=L, sigma=sigma)
    for pt in grafting_pts_ALPHA_P:
        ALPHA_PS.add(plates.add_tether(sticky_end='alpha_p', pos=pt))
    for pt in grafting_pts_BETA_P:
        BETA_PS.add(plates.add_tether(sticky_end='beta_p', pos=pt))

    # Set up its initial separation and energy scales (from here on,
    # the configurational binding entropies of all pairs of tethers
    # are fixed, so its much quicker to change binding energies)
    #
    # It doesn't matter what beta_DeltaG0 is set to at this stage, other
    # than it's something different from infinity for all the pairs
    # of strands that may bind during this run
    #
    print "Initializing configurational entropy factors..."
    plates.beta_DeltaG0['alpha', 'alpha_p'] = -10
    plates.beta_DeltaG0['alpha', 'beta_p'] = -5
    plates.beta_DeltaG0['beta', 'alpha_p'] = -5
    plates.separation = L
    plates.update()
    print "Done"


# For a given difference between strong and weak bond binding strength,
# map out bonding probabilities
def do_it(beta_DeltaDeltaG):

    print "Looking at beta_DeltaDeltaG = %g" % beta_DeltaDeltaG

    with open('competing_deltaDelta%g.txt' % beta_DeltaDeltaG, 'w') as f:
        for beta_DeltaG0Strong in xrange(0, -51, -1):
            print "    beta_DeltaG0Strong = %g" % beta_DeltaG0Strong
            beta_DeltaG0Weak = beta_DeltaG0Strong + beta_DeltaDeltaG

            plates.beta_DeltaG0['alpha', 'alpha_p'] = beta_DeltaG0Strong
            plates.beta_DeltaG0['alpha', 'beta_p'] = beta_DeltaG0Weak
            plates.beta_DeltaG0['beta', 'alpha_p'] = beta_DeltaG0Weak

            plates.update(DeltaG0_only=True)

            # Count strong and weak bonds
            num_strong_bonds = plates.count_bonds(ALPHAS, ALPHA_PS)
            num_weak_bonds = (plates.count_bonds(ALPHAS, BETA_PS) +
                              plates.count_bonds(BETAS, ALPHA_PS))

            # Output
            f.write('%.2f\t%.3f\t%.3f\n' %
                    (beta_DeltaG0Strong, num_strong_bonds, num_weak_bonds))

# Main module
setup()
do_it(beta_DeltaDeltaG=3)
do_it(beta_DeltaDeltaG=5)
do_it(beta_DeltaDeltaG=8)
do_it(beta_DeltaDeltaG=11)
do_it(beta_DeltaDeltaG=14)
