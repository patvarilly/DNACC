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

import dnacc
from dnacc import units, physics
from dnacc.units import nm
from dnacc.utils import pbc_delta

from ssDNA_stats import *

import numpy as np
from math import pi, sqrt, acos, cos, atan2, log, exp
import subprocess
import os
from contextlib import contextmanager

try:
    os.mkdir('results')
except Exception:
    pass

# ========================================================================
# Parameters
# ========================================================================

# Everything that Bortolo has calculated uses Kuhn lengths as the unit
# of length
l_Kuhn = 5 * nm

# Calculate a pair potential for Rogers & Crocker-style A and B spheres
R = 550 * nm

N_A_in_A = 4800.0
N_B_in_B = 4200.0
N_A_in_AB = 2400.0
N_B_in_AB = 3500.0

# From Dinamelt, nearest-neighbor rule, 125 mM NaCl
include_dangling_ends = False  # What we think Rogers and Crocker did
if include_dangling_ends:
    #
    # 5'-TAATGCCTGTCTACC-3'  with  5'-TGAGTTGCGGTAGAC-3'
    #            -------                      -------
    DeltaH0 = -56.7 * units.kcal / units.mol
    DeltaS0 = -166.0 * units.cal / units.mol / units.K
else:
    #
    # 5'-GTCTACC-3'  with  5'-GGTAGAC-3'
    #    -------              -------
    DeltaH0 = -47.8 * units.kcal / units.mol
    DeltaS0 = -139.6 * units.cal / units.mol / units.K

T_A_B_arr = (np.array((36.0, 35.0, 33.0, 32.0, 30.5)) + 273.15) * units.K
T_AB_AB_arr = (np.array((32.0, 29.0, 27.0, 25.0, 24.0)) + 273.15) * units.K


# ========================================================================
# Utility code
# ========================================================================
# Adapted from: http://software.clapper.org/grizzled-python/epydoc/ \
#  grizzled.os-pysrc.html#working_directory
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
# System setup
# ========================================================================

# Use explicit tethers in a plate of size boxL ** 2 to estimate plate-plate
# potential, then use Derjaguin approximation
boxL = 50 * l_Kuhn


def set_up_A_B():
    use_tethers_A_B_dat = True

    global boxL
    if use_tethers_A_B_dat:
        # Override the above to match tethers_A_B.dat dimensions
        boxL = 100 * l_Kuhn

    plates = dnacc.Plates(Lx=boxL, Ly=boxL, periodic=True,
                          tether_statistics=ExplicitSSDNAStatistics())

    if use_tethers_A_B_dat:
        with open(os.path.join('mc_results_bortolo', 'A_B_results',
                               'tethers_A_B.dat'), 'r') as f:
            for line in f:
                fields = line.split()
                pos = (float(fields[1]) * l_Kuhn, float(fields[2]) * l_Kuhn)
                if fields[0] == 'a':
                    plates.add_tether(pos=pos, sticky_end='A', plate='lower')
                else:
                    plates.add_tether(pos=pos, sticky_end='B', plate='upper')
    else:
        # Generate tether positions deterministically
        np.random.seed(1)

        N_A = int(N_A_in_A / (4 * pi * R ** 2) * boxL ** 2)
        N_B = int(N_B_in_B / (4 * pi * R ** 2) * boxL ** 2)

        for (x, y) in boxL * np.random.random_sample((N_A, 2)):
            plates.add_tether(pos=(x, y), sticky_end='A', plate='lower')
        for (x, y) in boxL * np.random.random_sample((N_B, 2)):
            plates.add_tether(pos=(x, y), sticky_end='B', plate='upper')

    plates.beta_DeltaG0['A', 'B'] = 0.0

    return plates


def set_up_A_B_mf():
    plates = dnacc.PlatesMeanField(
        tether_statistics=MeanFieldSSDNAStatistics())

    plates.add_tether_type(
        sigma=N_A_in_A / (4 * pi * R ** 2), sticky_end='A', plate='lower')
    plates.add_tether_type(
        sigma=N_B_in_B / (4 * pi * R ** 2), sticky_end='B', plate='upper')

    plates.beta_DeltaG0['A', 'B'] = 0.0

    return plates


def set_up_AB_AB():
    use_tethers_AB_AB_dat = True

    global boxL
    if use_tethers_AB_AB_dat:
        # Override the above to match tethers_A_B.dat dimensions
        boxL = 100 * l_Kuhn

    stats = ExplicitSSDNAStatistics()
    plates = dnacc.Plates(Lx=boxL, Ly=boxL,
                          tether_statistics=stats, periodic=True)

    if use_tethers_AB_AB_dat:
        with open(os.path.join('mc_results_bortolo', 'AB_AB_results',
                               'tethers_AB_AB.dat'), 'r') as f:

            # The file has four blocks of data; the first two specify the
            # grafting points of sticky_end='A' tethers on each plate,
            # the next two refer to sticky_end='B' tethers

            last_plate = None
            block_num = 0

            for line in f:
                fields = line.split()
                if fields[0] != last_plate:
                    last_plate = fields[0]
                    block_num += 1

                pos = (float(fields[1]) * l_Kuhn, float(fields[2]) * l_Kuhn)
                sticky_end = ('A' if block_num <= 2 else 'B')
                plate = ('lower' if fields[0] == 'a' else 'upper')
                plates.add_tether(pos=pos, sticky_end=sticky_end,
                                  plate=plate)

            # Print out stats
            print('%d A tethers on lower, %g expected' %
                  (len(plates.find_tethers(sticky_end='A', plate='lower')),
                   N_A_in_AB / (4 * pi * R ** 2) * boxL ** 2))
            print('%d B tethers on lower, %g expected' %
                  (len(plates.find_tethers(sticky_end='B', plate='lower')),
                   N_B_in_AB / (4 * pi * R ** 2) * boxL ** 2))
            print('%d A tethers on upper, %g expected' %
                  (len(plates.find_tethers(sticky_end='A', plate='upper')),
                   N_A_in_AB / (4 * pi * R ** 2) * boxL ** 2))
            print('%d B tethers on upper, %g expected' %
                  (len(plates.find_tethers(sticky_end='B', plate='upper')),
                   N_B_in_AB / (4 * pi * R ** 2) * boxL ** 2))

    else:
        # Generate tether positions deterministically
        np.random.seed(1)

        N_A = int(N_A_in_AB / (4 * pi * R ** 2) * boxL ** 2)
        N_B = int(N_B_in_AB / (4 * pi * R ** 2) * boxL ** 2)

        for (x, y) in boxL * np.random.random_sample((N_A, 2)):
            plates.add_tether(pos=(x, y), sticky_end='A', plate='lower')
        for (x, y) in boxL * np.random.random_sample((N_B, 2)):
            plates.add_tether(pos=(x, y), sticky_end='B', plate='lower')

        for (x, y) in boxL * np.random.random_sample((N_A, 2)):
            plates.add_tether(pos=(x, y), sticky_end='A', plate='upper')
        for (x, y) in boxL * np.random.random_sample((N_B, 2)):
            plates.add_tether(pos=(x, y), sticky_end='B', plate='upper')

    plates.beta_DeltaG0['A', 'B'] = 0.0

    return plates


def set_up_AB_AB_mf():
    plates = dnacc.PlatesMeanField(
        tether_statistics=MeanFieldSSDNAStatistics())

    plates.add_tether_type(
        sigma=N_A_in_AB / (4 * pi * R ** 2), sticky_end='A', plate='lower')
    plates.add_tether_type(
        sigma=N_B_in_AB / (4 * pi * R ** 2), sticky_end='B', plate='lower')
    plates.add_tether_type(
        sigma=N_A_in_AB / (4 * pi * R ** 2), sticky_end='A', plate='upper')
    plates.add_tether_type(
        sigma=N_B_in_AB / (4 * pi * R ** 2), sticky_end='B', plate='upper')

    plates.beta_DeltaG0['A', 'B'] = 0.0

    return plates


# ========================================================================
# Calculations
# ========================================================================

def adjust_potentials_ref(potentials, h_arr):
    """Set zero of potential to that found at the largest h."""
    max_h_index = np.argmax(h_arr)
    for V_arr in potentials.values():
        V_arr -= V_arr[max_h_index]
    return potentials


def calc_potentials(plates, dg0_list, h_arr):
    potentials = dict()
    for dg0 in dg0_list:
        potentials[dg0] = np.zeros(h_arr.shape)

    for (ih, h) in enumerate(h_arr):
        print("Calculating at h = %g nm" % (h / nm))
        plates.at(h)

        for dg0 in dg0_list:
            plates.beta_DeltaG0['A', 'B'] = dg0
            if isinstance(plates, dnacc.System):
                plates.update(DeltaG0_only=True)
                V = plates.free_energy / (plates.dims[0] * plates.dims[1])
            else:
                plates.update()
                V = plates.free_energy_density

            potentials[dg0][ih] = V

    # Set zero of potential to that found at the largest h
    max_h_index = np.argmax(h_arr)
    for V_arr in potentials.values():
        V_arr -= V_arr[max_h_index]
    return potentials


# Produce results
def do_potentials(name_x, name_y, plates, plates_mf, dg0_arr, h_arr):

    print("Calculating %s-%s potentials, explicit tethers..." %
          (name_x, name_y))
    potentials = calc_potentials(plates, dg0_arr, h_arr)

    print("Calculating %s-%s potentials, mean field..." %
          (name_x, name_y))
    potentials_mf = calc_potentials(plates_mf, dg0_arr, h_arr)

    with open(os.path.join(
            'results',
            '%s_%s_plate_potentials.txt' % (name_x, name_y)), 'w') as f:

        f.write('# DeltaG0 (kT)\t' 'h (nm)\t'
                'Explicit pair potential (kT / nm^2)\t'
                'Mean field pair potential (kT / nm^2)\n')

        for dg0 in dg0_arr:
            for (h, V, Vmf) in zip(h_arr, potentials[dg0],
                                   potentials_mf[dg0]):
                f.write("%g\t%g\t%g\t%g\n" % (dg0, h / nm, V, Vmf))
            f.write("\n")

    try:
        with working_directory('gp'):
            subprocess.call(['gnuplot', 'plot_%s_%s_plate_potentials.gp' %
                             (name_x, name_y)])
    except Exception:
        pass

    with open(os.path.join(
            'results',
            '%s_%s_sphere_potentials.txt' % (name_x, name_y)), 'w') as f:

        f.write('# DeltaG0 (kT)\t' 'h (nm)\t'
                'Explicit + Derjaguin Pair potential (kT)\t'
                'Mean field + Derjaguin Pair potential (kT)\n')

        for dg0 in dg0_arr:
            V_sphere = dnacc.calc_spheres_potential(
                h_arr, potentials[dg0], R)
            V_sphere_mf = dnacc.calc_spheres_potential(
                h_arr, potentials_mf[dg0], R)

            for (h, V, Vmf) in zip(h_arr, V_sphere, V_sphere_mf):
                f.write("%g\t%g\t%g\t%g\n" % (dg0, h / nm, V, Vmf))
            f.write("\n")

    try:
        with working_directory('gp'):
            subprocess.call(['gnuplot', 'plot_%s_%s_sphere_potentials.gp' %
                             (name_x, name_y)])
    except Exception:
        pass


def read_Qij_over_Qi_Qj(h, realisation):
    """Read Qij_over_Qi_Qj for a realisation of two spheres at a distance h.

    Qij_over_Qi_Qj / rho_0 = exp( - beta DeltaG^cnf_ij )
    """

    Qij_over_Qi_Qj = dnacc.utils.SymDict()
    with open(os.path.join(
            'mc_results_bortolo', 'A_B_results', 'spheres_mc_explicit_DGIJ',
            'hdist_%d.%d' % (int(h / l_Kuhn), realisation),
            'hyb1.dat'), 'r') as f:

        # First line tells you how many tethers there are on the left
        # and right spheres
        numA, numB = (int(_) for _ in f.next().split())

        # Then, there are numA blocks of the form:
        #  number of partners k
        #    partner 1     Q_ij / Q_i Q_j    error bar
        #    ...
        #    partner k     Q_ij / Q_i Q_j    error bar
        #
        # Remember, the indices are 1-based, not 0-based like we use
        # in Python
        for i in xrange(numA):
            num_partners = int(f.next())
            for jj in xrange(num_partners):
                fields = f.next().split()

                j = int(fields[0]) - 1 + numA
                Qij_over_Qi_Qj[i, j] = float(fields[1]) * (1.0 / l_Kuhn ** 3)

    return Qij_over_Qi_Qj


def read_Qi_over_Qinf(h, realisation):
    """Read Qi_over_Qinf for a realisation of two spheres at a distance h.

    Qi_over_Qinf = exp( - beta G^rep_i )
    """

    Qi_over_Qinf = list()

    with open(os.path.join(
            'mc_results_bortolo', 'A_B_results', 'spheres_mc_explicit_DGIJ',
            'hdist_%d.%d' % (int(h / l_Kuhn), realisation),
            'rep1.dat'), 'r') as f:

        numA = int(f.next())
        for i in xrange(numA):
            Qi_over_Qinf.append(float(f.next()))

    with open(os.path.join(
            'mc_results_bortolo', 'A_B_results', 'spheres_mc_explicit_DGIJ',
            'hdist_%d.%d' % (int(h / l_Kuhn), realisation),
            'rep2.dat'), 'r') as f:

        numB = int(f.next())
        for j in xrange(numB):
            Qi_over_Qinf.append(float(f.next()))

    return np.asarray(Qi_over_Qinf)


def do_explicit_spheres(prefix, dg0_arr):

    # Bortolo has prepared explicit DG_ij factors for explicit tethers
    # on spheres at various distances.  See what SCT predicts for numbers
    # of bonds and pair potentials

    results = {}

    realisation_arr = range(1, 16)
    # Skip the one realisation that didn't end in time
    realisation_arr.remove(11)

    h_arr = np.linspace(1 * l_Kuhn, 6 * l_Kuhn, 6)

    def calculate():
        # Only tethers on one sphere that may overlap the other sphere are
        # included in Bortolo's file.  For all the other tethers,
        # Qi_over_Qinf is a fixed number, which we take to be 0.199 from a
        # more detailed calculation elsewhere
        #
        # Getting this right is important to get a correct zero of free
        # energy.
        Nmax = 500
        Qi_over_Qinf_nonint = 0.1999

        for realisation in realisation_arr:
            print '-- realisation %d' % realisation
            results[realisation] = {}

            for h in h_arr:
                print 'h = %g nm' % (h / nm)

                results[realisation][h] = {}

                # Read
                Qij_over_Qi_Qj = read_Qij_over_Qi_Qj(h, realisation)
                Qi_over_Qinf = read_Qi_over_Qinf(h, realisation)
                N = len(Qi_over_Qinf)
                assert N <= Nmax

                # Now produce Boltzmann binding factors for various dg0
                for dg0 in dg0_arr:
                    boltz_bind = dict()
                    boltz_fac = exp(-dg0)
                    for key, value in Qij_over_Qi_Qj.iteritems():
                        boltz_bind[key] = boltz_fac * value / (1 * units.M)
                    boltz_bind_csr = dnacc.utils.csr_matrix_from_dict(
                        (N, N), boltz_bind)

                    # Set up DNACC system
                    spheres = dnacc.Generic(boltz_bind_csr)

                    # Gather results
                    binding_free_energy = spheres.binding_free_energy
                    rep_free_energy = (
                        -np.sum(np.log(Qi_over_Qinf)) -
                        (Nmax - N) * log(Qi_over_Qinf_nonint) +
                        Nmax * log(Qi_over_Qinf_nonint))

                    results[realisation][h][dg0] = dict(
                        avg_num_bonds=spheres.avg_num_bonds,
                        binding_free_energy=binding_free_energy,
                        rep_free_energy=rep_free_energy,
                        free_energy=binding_free_energy + rep_free_energy)

            # Set zero of free energy at last h
            #max_h = max(results[realisation].keys())
            #for dg0 in dg0_arr:
            #    these_results = results[realisation][max_h][dg0]
            #    bind_ref = these_results['binding_free_energy']
            #    rep_ref = these_results['rep_free_energy']
            #    ref = these_results['free_energy']
            #    for h in h_arr:
            #        these_results['free_energy'] -= ref
            #        these_results['binding_free_energy'] -= bind_ref
            #        these_results['rep_free_energy'] -= rep_ref

            # Print results out
            with open(os.path.join(
                    'results',
                    prefix + '_results_%d.txt' % realisation), 'w') as f:

                f.write('# dg0 (kT)\t' 'h (nm)\t' 'Total F (kT)\t'
                        'F_att (kT)\t' 'F_rep (kT)\t' 'avg_num_bonds\n')
                for dg0 in dg0_arr:
                    for h in h_arr:
                        R = results[realisation][h][dg0]
                        f.write('%g\t%g\t%g\t%g\t%g\t%g\n' %
                                (dg0, h / nm,
                                 R['free_energy'],
                                 R['binding_free_energy'],
                                 R['rep_free_energy'],
                                 R['avg_num_bonds'],
                                 ))
                    f.write('\n')

            # Reshuffle results for easier plotting
            with open(os.path.join(
                    'results',
                    prefix + '_results_%d_by_h.txt' % realisation), 'w') as f:

                f.write('# dg0 (kT)\t' 'h (nm)\t' 'Total F (kT)\t'
                        'F_att (kT)\t' 'F_rep (kT)\t' 'avg_num_bonds\n')
                for h in h_arr:
                    for dg0 in dg0_arr:
                        R = results[realisation][h][dg0]
                        f.write('%g\t%g\t%g\t%g\t%g\t%g\n' %
                                (dg0, h / nm,
                                 R['free_energy'],
                                 R['binding_free_energy'],
                                 R['rep_free_energy'],
                                 R['avg_num_bonds'],
                                 ))
                    f.write('\n')

    def read_in():
        for realisation in realisation_arr:
            results[realisation] = {}

            # Read results in
            with open(os.path.join(
                    'results',
                    prefix + '_results_%d.txt' % realisation), 'r') as f:

                # Skip first line
                f.next()

                for line in f:
                    fields = line.split()
                    if not fields:
                        continue

                    dg0 = float(fields[0])
                    h = float(fields[1]) * nm

                    h = h_arr[np.argmin(np.abs(np.array(h_arr) - h))]
                    dg0 = dg0_arr[np.argmin(np.abs(np.array(dg0_arr) - dg0))]

                    if h not in results[realisation]:
                        results[realisation][h] = {}
                    if dg0 not in results[realisation][h]:
                        results[realisation][h][dg0] = {}

                    R = results[realisation][h][dg0]
                    R['free_energy'] = float(fields[2])
                    R['binding_free_energy'] = float(fields[3])
                    R['rep_free_energy'] = float(fields[4])
                    R['avg_num_bonds'] = float(fields[5])

    def do_avg():
        # Average results over all realisations
        avg_results = {}

        for h in h_arr:
            avg_results[h] = {}
            for dg0 in dg0_arr:
                R = avg_results[h][dg0] = {}

                for name in ('free_energy', 'binding_free_energy',
                             'rep_free_energy', 'avg_num_bonds'):
                    R[name] = (
                        sum(results[realisation][h][dg0][name]
                            for realisation in realisation_arr) /
                        len(realisation_arr))

        with open(os.path.join('results',
                               prefix + '_avg_results.txt'), 'w') as f:
            f.write('# dg0 (kT)\t' 'h (nm)\t' 'Total F (kT)\t'
                    'F_att (kT)\t' 'F_rep (kT)\t' 'avg_num_bonds\n')
            for dg0 in dg0_arr:
                for h in h_arr:
                    R = avg_results[h][dg0]
                    f.write('%g\t%g\t%g\t%g\t%g\t%g\n' %
                            (dg0, h / nm,
                             R['free_energy'],
                             R['binding_free_energy'],
                             R['rep_free_energy'],
                             R['avg_num_bonds'],
                             ))
                f.write('\n')

        with open(os.path.join('results',
                               prefix + '_avg_results_by_h.txt'), 'w') as f:
            f.write('# dg0 (kT)\t' 'h (nm)\t' 'Total F (kT)\t'
                    'F_att (kT)\t' 'F_rep (kT)\t' 'avg_num_bonds\n')
            for h in h_arr:
                for dg0 in dg0_arr:
                    R = avg_results[h][dg0]
                    f.write('%g\t%g\t%g\t%g\t%g\t%g\n' %
                            (dg0, h / nm,
                             R['free_energy'],
                             R['binding_free_energy'],
                             R['rep_free_energy'],
                             R['avg_num_bonds'],
                             ))
                f.write('\n')

    # Main body
    calculate()
    read_in()
    do_avg()

    try:
        with working_directory('gp'):
            subprocess.call(['gnuplot', 'plot_%s_avg_num_bonds.gp' %
                             prefix])
    except Exception:
        pass

# ========================================================================
# MAIN PROGRAM
# ========================================================================

do_explicit_spheres(prefix='rogers_crocker_A_B',
                    dg0_arr=[(DeltaH0 - T * DeltaS0) / (physics.kB * T)
                             for T in T_A_B_arr])
do_explicit_spheres(prefix='general_A_B',
                    dg0_arr=np.linspace(-30.0, 0.0, 31))

do_potentials('A', 'B', set_up_A_B(), set_up_A_B_mf(),
              dg0_arr=[(DeltaH0 - T * DeltaS0) / (physics.kB * T)
                       for T in T_A_B_arr],
              h_arr=np.linspace(5 * nm, 50 * nm, 91)
              )

do_potentials('AB', 'AB', set_up_AB_AB(), set_up_AB_AB_mf(),
              dg0_arr=[(DeltaH0 - T * DeltaS0) / (physics.kB * T)
                       for T in T_AB_AB_arr],
              h_arr=np.linspace(5 * nm, 39.5 * nm, 70)
              )
