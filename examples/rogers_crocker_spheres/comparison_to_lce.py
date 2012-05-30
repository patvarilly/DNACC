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

import numpy as np
from math import pi, sqrt, acos, cos, atan2, log
import subprocess
import scipy.interpolate
import os
from contextlib import contextmanager

import dnacc
from dnacc.tether_statistics import RodsGraftedOnPlates
from dnacc.units import nm
from dnacc.utils import pbc_delta

from ssDNA_stats import *

try:
    os.mkdir('results')
except Exception:
    pass

# Everything that Bortolo has calculated uses Kuhn lengths as the unit
# of length
l_Kuhn = 5 * nm


# ========================================================================
# Utility code
# ========================================================================
# Adapted from: http://software.clapper.org/grizzled-python/epydoc/ \
#   grizzled.os-pysrc.html#working_directory
@contextmanager
def working_directory(directory):
    """This function is intended to be used as a ``with`` statement context
    manager. It allows you to replace code like this:

    .. python::

       original_directory = _os.getcwd()
       try:
           _os.chdir(some_dir)
           ... bunch of code ...
       finally:
           _os.chdir(original_directory)

    with something simpler:

    .. python ::

       from __future__ import with_statement
       from grizzled.os import working_directory

       with working_directory(some_dir):
           ... bunch of code ...

    :Parameters:
        directory : str
            directory in which to execute

    :return: yields the ``directory`` parameter
    """
    original_directory = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(original_directory)


# ========================================================================
# MAIN PROGRAM
# ========================================================================

# Read in tethers
tethers = np.loadtxt(os.path.join('mc_results_bortolo', 'sym_A_B_plates',
                                  'tethers.dat'),
                     dtype={'names': ('type', 'x', 'y'),
                            'formats': ('S1', 'd', 'd')},
                     skiprows=1)

# Set up plate system
boxL = 129.98 * l_Kuhn
for run in [dict(nfree=8, nbound=16),
            dict(nfree=8, nbound=15),
            ]:

    bridges_data_dir = 'DGbridge_%d_%d' % (run['nfree'], run['nbound'])
    outfile = os.path.join('results', 'results_%d_%d.txt' %
                           (run['nfree'], run['nbound']))
    dgcnffile = os.path.join('results',
                             ('dgcnf_%d_%d' % (run['nfree'], run['nbound']))
                             + '_h%.1f.txt')

    # Set up explicit tethers system
    stats = ExplicitSSDNAStatistics(bridges_data_dir=bridges_data_dir)
    plates = dnacc.Plates(Lx=boxL, Ly=boxL,
                          tether_statistics=stats, periodic=True)

    for t in tethers:
        if t['type'] == 'a':
            plates.add_tether(pos=(t['x'] * l_Kuhn, t['y'] * l_Kuhn),
                              sticky_end='alpha', plate='lower')
        else:
            plates.add_tether(pos=(t['x'] * l_Kuhn, t['y'] * l_Kuhn),
                              sticky_end='alphap', plate='upper')

    plates.beta_DeltaG0['alpha', 'alphap'] = -10

    # Set up mean field system
    meanfield_stats = MeanFieldSSDNAStatistics(stats)
    plates_meanfield = dnacc.PlatesMeanField(meanfield_stats)

    plates_meanfield.add_tether_type(sticky_end='alpha', plate='lower',
                                     sigma=len(tethers) / 2 / boxL ** 2)
    plates_meanfield.add_tether_type(sticky_end='alphap', plate='upper',
                                     sigma=len(tethers) / 2 / boxL ** 2)

    # Produce results
    with open(outfile, 'w') as f:
        f.write('# h (nm)\t' 'DeltaG0 (kT)\t'
                'avg_num_bonds (explicit)\t' 'avg_num_bonds (mean field)\n')

        #for h in (10 * nm, 15 * nm, 20 * nm, 25 * nm):
        for h in (15 * nm, 25 * nm):
            print outfile, h / nm
            plates.at(h)
            plates_meanfield.at(h)

            # Output DeltaG^cnf_ij matrix
            with open(dgcnffile % (h / l_Kuhn), 'w') as dg_f:
                dg_f.write('# i\t' 'j\t' 'beta * DeltaG^cnf_ij\n')
                for (i, j) in sorted(plates._boltz_binding_cnf.keys()):
                    dg_f.write('%d\t%d\t%.10g\n' %
                               (i, j,
                                -log(plates._boltz_binding_cnf[i, j])))

            for dg in np.linspace(-15.0, 0.0, 60 + 1):

                plates.beta_DeltaG0['alpha', 'alphap'] = dg
                plates.update(DeltaG0_only=True)

                plates_meanfield.beta_DeltaG0['alpha', 'alphap'] = dg
                plates_meanfield.update()

                f.write("%g\t%g\t%g\t%g\n" %
                        (h / nm, dg, plates.avg_num_bonds,
                         plates_meanfield.sigma_bonds * boxL ** 2))
                print("%g\t%g\t%g\t%g" %
                      (h / nm, dg, plates.avg_num_bonds,
                       plates_meanfield.sigma_bonds * boxL ** 2))

            f.write("\n")

    try:
        with working_directory('gp'):
            subprocess.call(['gnuplot', 'plot_%d_%d_lce_comparison.gp' %
                             (run['nfree'], run['nbound'])])
    except Exception:
        pass
