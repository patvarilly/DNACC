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

# Tether statistics (explicit and mean field) for ssDNA tethers on plates
# Written by Patrick Varilly, March 2012

# For a given ssDNA model (i.e., # segments when free, number of segments
# when bound), bridges_data_dir is a folder with a set of files named
# ('config_h%.1f.dat' % h), where h is the plate separation in units of Kuhn
# lengths (here, 5 nm, as in Rogers & Crocker, PNAS, 2011).  These files
# contain Q_ij / Q_i Q_j values for two tethers grafted between plates
# separated by h and where the vector joining the grafting points makes an
# angle theta with the normal of the plates.  From Bortolo's READ_ME.txt
# file, the format is:
#
# =================================================================
#
# FILE STRUCTURE (notice the first two lines)
#
# Number of tethered computed
# dtheta
#  ..................
# theta_i   ; Q_ij/Q_i Q_j (or exp(-\beta DGconf))
#  ..................
#
# theta_i goes from 0 to ACos[h/16] - dtheta
#
# Note: Q_ij has not been divided by rho_0 ( = 75.2768 * (5 nm)**3),
# which matches the definition of Q^cnf_ij in the paper.
#
# =================================================================
#
# For loops, the data is stored in loops_data_dir, and instead of
# parametrizing according to the angle theta, parametrize according to the
# distance s along the plate between the grafting points
#
# FILE STRUCTURE (notice the first two lines)
#
# Number of tethered computed
# ds
#  ..................
# s_i   ; Q_ij/Q_i Q_j (or exp(-\beta DGconf))
#  ..................
#
# s_i goes from 0 + ds/2 to 16 - ds/2
#
# =================================================================
#
# Finally, rep_data_file contains a table of entries (h, Qh_over_Qinf) that
# quantify the entropic cost of confining a free tether between two plates
# at a distance h
#
# Note: Q_ij has not been divided by rho_0 ( = 75.2768 * (5 nm)**3),
# which matches the definition of Q^cnf_ij in the paper.

__all__ = ['ExplicitSSDNAStatistics', 'MeanFieldSSDNAStatistics']

import os

data_dir = os.path.join('mc_results_bortolo', 'ssDNA_plates_DGIJ')

import numpy as np
from math import pi, sqrt, acos, cos, atan2, log, exp
import scipy.interpolate

import dnacc
from dnacc.tether_statistics import RodsGraftedOnPlates
from dnacc import units
from dnacc.units import nm
from dnacc.utils import pbc_delta

# This example needs SciPy 0.9.0 or later for robust 2D interpolation
# on an irregular grid
import scipy
majver, minver = scipy.__version__.split('.')[:2]
assert majver >= 0 or minver >= 9


class ExplicitSSDNAStatistics(object):
    def __init__(self,
                 bridges_data_dir='DGbridge_8_16',
                 loops_data_dir='DGloop_8_16',
                 rep_data_file='QhOverQinf.dat'):

        bridges_data_dir = os.path.join(data_dir, bridges_data_dir)
        loops_data_dir = os.path.join(data_dir, loops_data_dir)
        rep_data_file = os.path.join(data_dir, rep_data_file)

        # Everything that Bortolo has calculated uses Kuhn lengths as the
        # unit of length.  We convert to DNACC units as soon as we read
        # data in
        l_Kuhn = 5 * nm

        # Repulsion
        self.Qh_over_Qinf_data = np.loadtxt(rep_data_file)
        self.Qh_over_Qinf_data[:, 0] *= l_Kuhn

        self.Qh_over_Qinf = scipy.interpolate.interp1d(
            self.Qh_over_Qinf_data[:, 0], self.Qh_over_Qinf_data[:, 1],
            kind='cubic', bounds_error=True)

        # Bridges
        self.bridge_points = []
        self.bridge_data = []
        for h in np.linspace(1.0, 10.0, 91):
            data_h = np.loadtxt(os.path.join(bridges_data_dir,
                                             'config_h%.1f.dat' % h),
                                skiprows=2)

            # Points are sampled at i * dtheta

            self.bridge_points.extend((h * l_Kuhn, theta)
                                      for theta in data_h[:, 0])
            self.bridge_data.extend(data_h[:, 1] * (1.0 / l_Kuhn ** 3))

            # Add final data point
            dtheta = data_h[1, 0] - data_h[0, 0]
            self.bridge_points.append((h * l_Kuhn, data_h[-1, 0] + dtheta))
            self.bridge_data.append(0.0)

        self.bridge_points = np.array(self.bridge_points)
        self.bridge_data = np.array(self.bridge_data)

        self.Qij_over_Qi_Qj_bridge = (
            scipy.interpolate.CloughTocher2DInterpolator(
                self.bridge_points, self.bridge_data, fill_value=0.0))

        # Loops
        if not os.path.exists(loops_data_dir):
            self.loop_points = None
            self.loop_data = None
            self.Qij_over_Qi_Qj_loop = None
        else:
            self.loop_points = []
            self.loop_data = []
            for h in np.linspace(0.2, 7.9, 78):
                data_h = np.loadtxt(os.path.join(loops_data_dir,
                                                 'dist_%g' % h,
                                                 'config.dat'),
                                    skiprows=2)

                # Points are sampled at (i+0.5) * ds

                # Add initial data point
                ds = (data_h[1, 0] - data_h[0, 0]) * l_Kuhn
                self.loop_points.append((h * l_Kuhn, 0.0))
                self.loop_data.append(
                    (data_h[0, 1] - 0.5 * (data_h[1, 1] - data_h[0, 1])) *
                    (1.0 / l_Kuhn ** 3))

                # Add bulk of data
                self.loop_points.extend((h * l_Kuhn, s * l_Kuhn)
                                        for s in data_h[:, 0])
                self.loop_data.extend(data_h[:, 1] * (1.0 / l_Kuhn ** 3))

                # Add final data point
                self.loop_points.append((h * l_Kuhn,
                                         data_h[-1, 0] * l_Kuhn + 0.5 * ds))
                self.loop_data.append(0.0)

            self.loop_points = np.array(self.loop_points)
            self.loop_data = np.array(self.loop_data)

            self.Qij_over_Qi_Qj_loop = (
                scipy.interpolate.CloughTocher2DInterpolator(
                    self.loop_points, self.loop_data, fill_value=0.0))

        # Find largest distance between the grafting points of two tethers
        # that can interact
        self.min_h_rep = min(self.Qh_over_Qinf_data[:, 0])
        self.max_h_rep = max(self.Qh_over_Qinf_data[:, 0])

        self.min_h_bridge = min(
            h for (h, theta) in self.bridge_points)
        self.max_h_bridge = max(
            h for (h, theta) in self.bridge_points)
        self.max_theta_bridge = max(
            theta for (h, theta) in self.bridge_points)
        self.max_d_bridge = max(
            h / cos(theta) for (h, theta) in self.bridge_points)

        if self.Qij_over_Qi_Qj_loop:
            self.min_h_loop = min(h for (h, s) in self.loop_points)
            self.max_h_loop = max(h for (h, s) in self.loop_points)
            self.max_s_loop = max(s for (h, s) in self.loop_points)
            self.max_d_loop = self.max_s_loop
        else:
            self.max_d_loop = 0.0

        self.max_d = max((self.max_d_bridge, self.max_d_loop))

    def calc_interaction_lists(self, system):
        return (RodsGraftedOnPlates().
                _cell_list_interaction_lists(system, self.max_d))

    def calc_boltz_binding_cnf(self, system, tether_info_i, tether_info_j):
        h = system.separation
        ri = tether_info_i['pos']
        rj = tether_info_j['pos']

        if system.periodic:
            boxL = system.dims
            dr = (pbc_delta(ri[0], rj[0], boxL[0]),
                  pbc_delta(ri[1], rj[1], boxL[1]))
        else:
            dr = (ri[0] - rj[0], ri[1] - rj[1])

        if tether_info_i['plate'] != tether_info_j['plate']:
            # Calculate angle between vector joining grafting points and
            # plate normal
            theta = atan2(sqrt(dr[0] ** 2 + dr[1] ** 2), h)
            Qij_over_Qi_Qj = max((0.0,
                                  self.Qij_over_Qi_Qj_bridge(h, theta)))

        else:
            if not self.Qij_over_Qi_Qj_loop:
                raise NotImplementedError('No loop data')

            # Watch out! Data goes out to h = self.max_h_loop, and then
            # the interpolator returns Qij_over_Qi_Qj == 0.0, which is not
            # true.  Guard explicitly against that possibility
            h = min(0.99 * self.max_h_loop, h)
            s = sqrt(dr[0] ** 2 + dr[1] ** 2)
            Qij_over_Qi_Qj = max((0.0, self.Qij_over_Qi_Qj_loop(h, s)))

        Qi_over_Qinf = self.calc_boltz_exclusion(system, tether_info_i)
        # Optimization
        Qj_over_Qinf = Qi_over_Qinf

        return (Qij_over_Qi_Qj * Qi_over_Qinf * Qj_over_Qinf /
                (1 * units.M))

    def calc_boltz_exclusion(self, system, tether_info_i):
        h = system.separation
        if h > self.max_h_rep:
            return 1.0
        elif h < self.min_h_rep:
            return 0.0
        else:
            return max(0.0, self.Qh_over_Qinf(h))


class MeanFieldSSDNAStatistics(object):
    # We compute a mean field K_{ij} by direct integration of the
    # Boltzmann factors exp(-beta Delta G^cnf_{ij})

    def __init__(self, explicit=None):
        if explicit:
            self._explicit = explicit
        else:
            self._explicit = ExplicitSSDNAStatistics()

        # We cache results according to plate separation
        self._bridge_boltz = dict()
        self._loop_boltz = dict()

    def _calc(self, h, kind):
        if not (kind == 'bridge' or kind == 'loop'):
            raise ValueError('Unrecognized kind %s' % str(kind))

        if kind == 'bridge':
            min_d, max_d = h, self._explicit.max_d_bridge
        else:
            min_d, max_d = 0.0, self._explicit.max_d_loop

        if max_d < min_d:
            return 0.0

        class FakeSystem(object):
            pass
        system = FakeSystem()
        system.separation = h
        system.periodic = False

        if kind == 'bridge':
            def integrand(d):
                s = max(0.0, sqrt(d ** 2 - h ** 2))
                return d * self._explicit.calc_boltz_binding_cnf(
                    system,
                    dict(pos=(0.0, 0.0), plate='lower'),
                    dict(pos=(0.0, s), plate='upper'))
        else:
            def integrand(d):
                return d * self._explicit.calc_boltz_binding_cnf(
                    system,
                    dict(pos=(0.0, 0.0), plate='lower'),
                    dict(pos=(0.0, d), plate='lower'))

        epsabs = 1e-5 * nm ** 2
        epsrel = 1e-5

        return (2 * pi * scipy.integrate.quad(integrand, min_d, max_d,
                                              epsabs=epsabs,
                                              epsrel=epsrel)[0])

    def calc_boltz_exclusion(self, system, type_info_a):
        return self._explicit.calc_boltz_exclusion(system, type_info_a)

    def calc_boltz_binding_cnf_bridge(self, system, type_info_a, type_info_b):
        h = system.separation
        if h not in self._bridge_boltz:
            value = self._bridge_boltz[h] = self._calc(h, kind='bridge')
        else:
            value = self._bridge_boltz[h]
        return value

    def calc_boltz_binding_cnf_loop(self, system, type_info_a, type_info_b):
        h = system.separation
        if h not in self._loop_boltz:
            value = self._loop_boltz[h] = self._calc(h, kind='loop')
        else:
            value = self._loop_boltz[h]
        return value

# Testing
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    if False:

        explicit = ExplicitSSDNAStatistics()

        #fig = plt.figure()
        #ax = Axes3D(fig)
        #ax.scatter3D([h / nm for (h, theta) in explicit.bridge_points],
        #             [theta for (h, theta) in explicit.bridge_points],
        #             explicit.bridge_data)
        #plt.show()

        #grid_x, grid_y = np.mgrid[0 : explicit.max_h_bridge     : 100j,
        #                          0 : explicit.max_theta_bridge : 100j]
        #fig = plt.figure()
        #ax = Axes3D(fig)
        #ax.scatter3D(grid_x.ravel(), grid_y.ravel(),
        #             explicit.Qij_over_Qi_Qj_bridge((grid_x, grid_y)).ravel())
        #plt.show()

        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter3D([h / nm for (h, s) in explicit.loop_points],
                     [s / nm for (h, s) in explicit.loop_points],
                     explicit.loop_data)
        plt.show()

        grid_x, grid_y = np.mgrid[0:explicit.max_h_loop:100j,
                                  0:explicit.max_s_loop:100j]
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter3D(grid_x.ravel(), grid_y.ravel(),
                     explicit.Qij_over_Qi_Qj_loop((grid_x, grid_y)).ravel())
        plt.show()

    elif False:

        mf = MeanFieldSSDNAStatistics()

        h_arr = np.linspace(4 * nm, mf._explicit.max_h_rep, 50)

        print 'Calculating K^bridge(h)...'
        fig = plt.figure()
        plt.plot(h_arr, [mf._calc(h, kind='bridge') / nm ** 2
                         for h in h_arr])
        plt.show()

        print 'Calculating K^loop(h)...'
        fig = plt.figure()
        plt.plot(h_arr, [mf._calc(h, kind='loop') / nm ** 2
                         for h in h_arr])
        plt.show()

        print 'Done'

    else:

        explicit = ExplicitSSDNAStatistics()

        l_Kuhn = 5 * nm
        h = 0.5 * l_Kuhn
        s_arr = np.linspace(0 * nm, 16 * l_Kuhn, 200)
        Qij_over_Qi_Qj = [explicit.Qij_over_Qi_Qj_loop((h, s))
                          for s in s_arr]
        plt.plot(s_arr, Qij_over_Qi_Qj)

        with open('temp.txt', 'w') as f:
            for (s, Q) in zip(s_arr, Qij_over_Qi_Qj):
                f.write('%g\t%g\n' % (s, Q))

        plt.show()
