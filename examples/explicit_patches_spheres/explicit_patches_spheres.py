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


def setup_explicit(sphere_list, patches):

    # Set up system
    spheres = dnacc.Spheres()
    spheres.set_tether_prototype(L=L)

    for name, props in sphere_list.iteritems():
        spheres.add_sphere(name, props['centre'], props['radius'])

    # Generate grafting points for patches centered at 0.0
    ids = list()
    for patch_props in patches:

        patch_ids = dnacc.patches.add_circular_patch_to_sphere(
            spheres, **patch_props)

        ids.extend(patch_ids)

    # Record original positions
    for i in ids:
        spheres.tethers[i]['orig_pos'] = np.array(spheres.tethers[i]['pos'])

    return spheres


# Main module
# ===========

np.random.seed(1)  # To make the following plots deterministic

L = 20 * nm

system = setup_explicit(

    sphere_list=dict(
        left_ball=dict(centre=(0.0, 0.0, 0.0), radius=500 * nm),
        right_ball=dict(centre=(1000.0 * nm, 0.0, 0.0), radius=500 * nm)),

    patches=[
        dict(sphere='left_ball', centre=(500.0 * nm, 0.0, 0.0),
             angular_aperture=pi / 24.0, N=10, L=L,
             sticky_end='alpha'),

        dict(sphere='right_ball', centre=(-500.0 * nm, 0.0, 0.0),
             angular_aperture=pi / 24.0, N=10, L=L,
             sticky_end='alphap')])

system.beta_DeltaG0['alpha', 'alphap'] = -15

with open('results.txt', 'w') as f:
    print('# Separation (nm)\t' 'F (kT)')
    f.write('# Separation (nm)\t' 'F (kT)\n')
    for h in np.linspace(0.05 * L, 2.00 * L, 40):
        system.sphere_centre_at('right_ball', (1000 * nm + h, 0.0, 0.0))
        F = system.free_energy
        print "%g\t%g" % (h / nm, F)
        f.write("%g\t%g\n" % (h / nm, F))

subprocess.call(['gnuplot', 'plot_results.gp'])
