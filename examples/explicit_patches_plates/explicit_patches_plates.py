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
# Python script to generate interaction potentials for two patches
# on opposing parallel plates

import numpy as np
from math import pi, sqrt
import subprocess
import scipy.interpolate

import dnacc
from dnacc.units import nm


def setup_explicit(L, sigma, patch_radii):
    # Set up system
    boxL = 4 * max(patch_radii.itervalues()) + 4 * L
    plates = dnacc.Plates(boxL, boxL, periodic=True)
    plates.set_tether_prototype(L=L)

    # Generate grafting points for patches centered at 0.0
    ids = list()
    for patch_props in (dict(plate='upper', sticky_end='alpha'),
                        dict(plate='lower', sticky_end='alphap')):

        R = patch_radii[patch_props['plate']]
        N = int(pi * R ** 2 * sigma)
        patch_ids = dnacc.patches.add_circular_patch_to_plate(
            plates, centre=(0.0, 0.0), R=R, N=N, **patch_props)

        ids.extend(patch_ids)

    # Record original positions
    for i in ids:
        plates.tethers[i]['orig_pos'] = plates.tethers[i]['pos']

    return plates


def calc_explicit_potential(plates, beta_DeltaG0, h, d):
    plates.separation = h
    plates.beta_DeltaG0['alpha', 'alphap'] = beta_DeltaG0

    for tether in plates.tethers:
        if tether['plate'] == 'upper':
            tether['pos'] = tether['orig_pos'] + np.array((d, 0.0))

    plates.update()

    return plates.free_energy


def run(L, sigma, beta_DeltaG0, patch_radii, skip=False, plates=None):

    if plates is None:
        plates = setup_explicit(L, sigma, patch_radii)
    max_patch_radii = max(patch_radii.itervalues())

    if skip:
        return plates

    filename = ("patch_results_L%gnm_S%g_dg%gkT_upperR%gnm_lowerR%gnm.dat" %
                (L / nm, sqrt(1 / sigma) / L, beta_DeltaG0,
                 patch_radii['upper'] / nm, patch_radii['lower'] / nm))
    print "Working on %s" % filename
    print "Total number of tethers: %d" % len(plates.tethers)

    with open(filename, 'w') as f:
        f.write("# d (nm)\t" "h (nm)\t" "F (kT)\n")
        f.write("# Total number of tethers: %d\n" % len(plates.tethers))
        for d in np.linspace(-2 * max_patch_radii - 2 * L,
                             +2 * max_patch_radii + 2 * L, 40):
            for h in np.linspace(0.05 * L, 2 * L, 20):
                f.write("%g\t%g\t%g\n" %
                        (d / nm, h / nm,
                         calc_explicit_potential(plates, beta_DeltaG0,
                                                 h, d)))

            f.write("\n")

    return plates

# Main module
# ===========

np.random.seed(1)  # To make the following plots deterministic

run(L=20 * nm, sigma=1 / (20 * nm) ** 2, beta_DeltaG0=-10,
    patch_radii=dict(upper=100 * nm, lower=100 * nm), skip=True)
subprocess.call(['gnuplot', 'plot_patch1.gp'])

run(L=20 * nm, sigma=1 / (20 * nm) ** 2, beta_DeltaG0=-10,
    patch_radii=dict(upper=100 * nm, lower=50 * nm), skip=True)
subprocess.call(['gnuplot', 'plot_patch2.gp'])

run(L=10 * nm, sigma=1 / (10 * nm) ** 2, beta_DeltaG0=-10,
    patch_radii=dict(upper=10 * nm, lower=10 * nm), skip=True)
subprocess.call(['gnuplot', 'plot_patch3.gp'])

run(L=10 * nm, sigma=1 / (8 * nm) ** 2, beta_DeltaG0=-10,
    patch_radii=dict(upper=10 * nm, lower=10 * nm), skip=True)
subprocess.call(['gnuplot', 'plot_patch4.gp'])

run(L=10 * nm, sigma=1 / (5 * nm) ** 2, beta_DeltaG0=-10,
    patch_radii=dict(upper=10 * nm, lower=10 * nm), skip=True)
subprocess.call(['gnuplot', 'plot_patch5.gp'])

nice_patch = run(L=20 * nm, sigma=1 / (20 * nm) ** 2, beta_DeltaG0=-8,
                 patch_radii=dict(upper=100 * nm, lower=100 * nm))
run(L=20 * nm, sigma=1 / (20 * nm) ** 2, beta_DeltaG0=-7,
    patch_radii=dict(upper=100 * nm, lower=100 * nm),
    plates=nice_patch)
run(L=20 * nm, sigma=1 / (20 * nm) ** 2, beta_DeltaG0=-6,
    patch_radii=dict(upper=100 * nm, lower=100 * nm),
    plates=nice_patch)
subprocess.call(['gnuplot', 'plot_patch6.gp'])
