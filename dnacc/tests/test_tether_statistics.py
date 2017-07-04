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

import unittest
import numpy as np
from dnacc.tether_statistics import *
from dnacc.units import nm
from math import *
import scipy

from nose.tools import *


## Unit testing
## ------------

# Rods grafted on plates
class TestRodsGraftedOnPlates(unittest.TestCase):
    """Unit test for RodsGraftedOnPlates."""

    def testSymmetricCase(self):
        calc = RodsGraftedOnPlates()

        La = Lb = 20 * nm

        # Test when vector joining anchors is perpendicular to plates
        for h in np.linspace(-0.5 * nm, La + Lb + 0.5 * nm, 50):
            if h <= 0.0:
                expOmegaA = expOmegaB = 0.0
            else:
                expOmegaA = 2 * pi * La * min(La, h)
                expOmegaB = 2 * pi * Lb * min(Lb, h)

            if not 0.0 < h < La + Lb:
                expOmegaAB = 0.0
            else:
                gammaA = acos((La ** 2 - Lb ** 2 + h ** 2) /
                              (2 * h * La))
                expOmegaAB = (2 * pi * La * sin(gammaA) /
                              sin(2 * gammaA))

            self.assertAlmostEqual(calc._Omega_bridge(La, Lb, h, h),
                                   expOmegaAB)
            self.assertAlmostEqual(calc._Omega_unbound(La, h),
                                   expOmegaA)
            self.assertAlmostEqual(calc._Omega_unbound(Lb, h),
                                   expOmegaB)

    def test_a_few_tethers(self):
        # Compare a set of real-world tethers to what
        # comes out of my older C++ code, which is thoroughly debugged
        calc = RodsGraftedOnPlates()

        test_cases = [
            # (Li, Lj, h, d, result)  <-- units are nm
            (20, 20, 2, 18.352, 0.000129815),
            (20, 20, 2, 18.7676, 0.00012769),
            (20, 20, 2, 28.6733, 0.000105543),
            (20, 20, 2, 12.191, 0.000183705),
            (20, 20, 2, 12.7734, 0.000175997),
            (20, 20, 2, 33.4347, 0.000114955),
            (20, 20, 2, 13.5029, 0.000167391),
            (20, 20, 2, 29.5537, 0.000105941),
            (20, 20, 2, 14.3053, 0.000159066),
            (20, 20, 2, 27.1822, 0.000105832),
            (20, 20, 2, 36.7817, 0.000146106),
            (20, 20, 2, 31.9561, 0.000109763),
            (20, 20, 2, 14.8298, 0.000154181),
            (20, 20, 2, 11.0637, 0.000201215),
            (20, 20, 2, 9.44597, 0.000234549),
            (20, 20, 2, 17.8691, 0.000132453),
            (20, 20, 2, 29.2586, 0.000105757),
            (20, 20, 2, 21.9646, 0.000115115),
            (20, 20, 2, 33.1272, 0.000113632),
            (20, 20, 2, 35.5812, 0.000129838),
            (20, 20, 2, 29.2128, 0.000105733),
            (20, 20, 2, 34.8743, 0.000123551),
            (20, 20, 2, 18.3633, 0.000129756),
            (20, 20, 2, 38.0214, 0.000179102),
            (20, 20, 2, 34.5325, 0.000121079),
            (20, 20, 2, 29.5489, 0.000105938),
            (20, 20, 2, 28.2874, 0.000105507),
            (20, 20, 2, 32.1512, 0.000110295),
            (20, 20, 2, 17.1877, 0.000136506),
            (20, 20, 2, 36.3728, 0.000139508),
            (20, 20, 2, 15.5233, 0.000148311),
            (20, 20, 2, 35.8734, 0.000133019),
            (20, 20, 2, 32.9492, 0.000112931),
            (20, 20, 2, 37.1601, 0.000153624),
            (20, 20, 2, 21.5643, 0.000116365),
            (20, 20, 2, 33.8003, 0.00011673),
            (20, 20, 2, 36.1966, 0.000137046),
            (20, 20, 2, 27.9447, 0.000105542),
            (20, 20, 2, 33.3591, 0.000114616),
            (20, 20, 2, 6.51025, 0.00034421),
            (20, 20, 2, 12.7578, 0.000176193),
            (20, 20, 2, 15.3329, 0.00014986),
            (20, 20, 2, 13.1086, 0.000171908),
            (20, 20, 2, 21.3819, 0.000116963),
            (20, 20, 2, 19.9781, 0.000122193),
            (20, 20, 2, 22.5105, 0.000113542),
            (20, 20, 2, 35.3381, 0.000127472),
            (20, 20, 2, 26.4642, 0.000106364),
            (20, 20, 2, 21.0651, 0.000118045),
            (20, 20, 2, 25.7753, 0.000107099),
            (20, 20, 2, 36.1867, 0.000136914),
            (20, 20, 2, 34.7924, 0.000122929),
            (20, 20, 2, 31.1747, 0.000108006),
            (20, 20, 2, 17.2531, 0.000136099),
            (20, 20, 2, 27.1591, 0.000105846),
            (20, 20, 2, 33.2938, 0.00011433),
            (20, 20, 2, 27.6373, 0.000105624),
            (20, 20, 2, 32.739, 0.000112159),
            (20, 20, 2, 23.5751, 0.000110897),
            (20, 20, 2, 23.8693, 0.000110261),
            (20, 20, 2, 34.9251, 0.000123947),
            (20, 20, 2, 20.3142, 0.000120835),
            (20, 20, 2, 25.4674, 0.000107498),
            (20, 20, 2, 25.301, 0.000107731),
            (20, 20, 2, 25.8851, 0.000106967),
            (20, 20, 2, 35.6605, 0.000130662),
            (20, 20, 2, 38.22, 0.000187699),
            (20, 20, 2, 26.7343, 0.000106135),
            (20, 20, 2, 36.2144, 0.000137286),
            (20, 20, 2, 34.0973, 0.000118355),
            (20, 20, 2, 29.2384, 0.000105746),
            (20, 20, 2, 36.0173, 0.00013474),
            (20, 20, 2, 34.5037, 0.000120885),
            (20, 20, 2, 20.1534, 0.000121477),
            (20, 20, 2, 23.2669, 0.000111607),
            (20, 20, 2, 38.4582, 0.000200267),
            (20, 20, 2, 34.2417, 0.000119211),
            (20, 20, 2, 32.8076, 0.000112405),
            (20, 20, 2, 9.19438, 0.000240908),
            (20, 20, 2, 17.8882, 0.000132345),
            (20, 20, 2, 31.4721, 0.000108608),
            (20, 20, 2, 33.1687, 0.000113802)]

        for case in test_cases:
            self.assertAlmostEqual(
                calc._calc_bridge_boltz(
                    case[0] * nm, case[1] * nm,
                    case[2] * nm, case[3] * nm) /
                (calc._calc_boltz_exclusion(case[0] * nm, case[2] * nm) *
                 calc._calc_boltz_exclusion(case[1] * nm, case[2] * nm)),
                case[4])

    def test_distance(self):
        # Unit test to capture bug in calculating distance between
        # two tethering points on different plates
        class FakeSystem(object):
            pass

        system = FakeSystem()
        system.dims = (200 * nm, 200 * nm)
        system.periodic = True

        tether_on_A = dict(pos=(0.0, 0.0), L=20 * nm, plate='A')
        tether_on_B = dict(pos=(sqrt((18.352 * nm) ** 2 -
                                     (2 * nm) ** 2), 0.0),
                           L=20 * nm, plate='B')
        tether2_on_A = dict(pos=(20 * nm, 0.0), L=20 * nm, plate='A')
        system.separation = 2 * nm

        stats = RodsGraftedOnPlates()
        stats.check_system(system)  # No exception here!

        self.assertAlmostEqual(
            stats.calc_boltz_binding_cnf(system, tether_on_A, tether_on_B) /
            (stats.calc_boltz_exclusion(system, tether_on_A) *
             stats.calc_boltz_exclusion(system, tether_on_B)),
            0.000129815)

        self.assertAlmostEqual(
            stats.calc_boltz_binding_cnf(system, tether_on_A, tether2_on_A) /
            (stats.calc_boltz_exclusion(system, tether_on_A) *
             stats.calc_boltz_exclusion(system, tether2_on_A)),
            0.00012169381364824528)

    def test_checks(self):
        class FakeSystem(object):
            pass

        stats = RodsGraftedOnPlates()

        # Partially defined systems
        system = FakeSystem()
        self.assertRaises(ValueError, stats.check_system, system)

        system.dims = (20 * nm, 20 * nm)
        self.assertRaises(ValueError, stats.check_system, system)

        system.periodic = True
        self.assertRaises(ValueError, stats.check_system, system)

        system.separation = 10 * nm
        stats.check_system(system)  # No exception here!

        # Invalid inputs
        system.dims = (20 * nm)
        self.assertRaises(ValueError, stats.check_system, system)
        system.dims = (20 * nm, 20 * nm, 20 * nm)
        self.assertRaises(ValueError, stats.check_system, system)
        system.dims = (20 * nm, 20 * nm)

        system.separation = 0
        self.assertRaises(ValueError, stats.check_system, system)
        system.separation = -10 * nm
        self.assertRaises(ValueError, stats.check_system, system)
        system.separation = 10 * nm

        stats.check_system(system)  # No exception here!

        # Now run checks on tethers

        tether = {}
        self.assertRaises(ValueError, stats.check_tether,
                          system, tether)

        tether['plate'] = 'upper'
        self.assertRaises(ValueError, stats.check_tether,
                          system, tether)

        tether['L'] = 20 * nm
        self.assertRaises(ValueError, stats.check_tether,
                          system, tether)

        tether['pos'] = 'hello'
        self.assertRaises(ValueError, stats.check_tether,
                          system, tether)

        tether['pos'] = (10 * nm,)
        self.assertRaises(ValueError, stats.check_tether,
                          system, tether)

        tether['pos'] = (10 * nm, 20 * nm, 30 * nm)
        self.assertRaises(ValueError, stats.check_tether,
                          system, tether)

        tether['pos'] = (10 * nm, 20 * nm)
        stats.check_tether(system, tether)  # No exception here!

    # Test interaction lists
    def test_interaction_lists(self):
        class FakeSystem(object):
            pass

        L = 20 * nm
        S = 0.75

        sigma = 1 / (S * L) ** 2

        N = 100

        boxL = sqrt(N / sigma)

        system = FakeSystem()
        system.dims = (boxL, boxL)
        system.periodic = True
        system.separation = 10 * nm

        np.random.seed(0)
        system.tethers = (
            list(
                dict(pos=(x * boxL, y * boxL), L=L, plate='lower')
                for x, y in np.random.random_sample((N, 2))) +
            list(
                dict(pos=(x * boxL, y * boxL), L=L, plate='upper')
                for x, y in np.random.random_sample((N, 2))))

        stats = RodsGraftedOnPlates()
        partners = stats.calc_interaction_lists(system)

        # Check that all non-partners far away
        boxL = system.dims
        if system.periodic:
            def distance(ri, rj):
                dx = ri[0] - rj[0]
                dx -= round(dx / boxL[0]) * boxL[0]
                dy = ri[1] - rj[1]
                dy -= round(dy / boxL[1]) * boxL[1]
                return sqrt(dx ** 2 + dy ** 2)
        else:
            def distance(ri, rj):
                return sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2)

        all_tethers = set(range(N))
        #for i, nns_i in enumerate(partners):
        #    for j in nns_i:
        #        assert (distance(system.tethers[i]['pos'],
        #                         system.tethers[j]['pos']) < 2 * L)

        for i, nns_i in enumerate(partners):
            for j in (all_tethers - set(nns_i)):
                assert (distance(system.tethers[i]['pos'],
                                 system.tethers[j]['pos']) >= 2 * L)


# Rods grafted on plates, mean field
class TestRodsGraftedOnPlatesMeanField(unittest.TestCase):
    """Unit test for RodsGraftedOnPlatesMeanField."""

    def test_equivalence(self):
        stats = RodsGraftedOnPlates()
        mean_field_stats = RodsGraftedOnPlatesMeanField()

        def test_it(La, Lb, h, loop=False, bridge=False):
            assert bool(bridge) ^ bool(loop)

            if bridge:
                minD = h
                maxD = La + Lb
            else:
                minD = abs(La - Lb)
                maxD = La + Lb

            if not (minD < maxD):
                return

            if bridge:
                def integrand(d):
                    return d * stats._calc_bridge_boltz(La, Lb, h, d)
            else:
                def integrand(d):
                    return d * stats._calc_loop_boltz(La, Lb, h, d)

            epsabs = 1e-5 * nm ** 2
            epsrel = 1e-5

            explicit = (
                2 * pi * scipy.integrate.quad(
                    integrand, minD, maxD,
                    epsabs=epsabs, epsrel=epsrel)[0] /
                (stats._calc_boltz_exclusion(La, h) *
                 stats._calc_boltz_exclusion(Lb, h)))

            class FakeSystem(object):
                pass
            system = FakeSystem()
            system.separation = h
            tether_type_a = dict(L=La)
            tether_type_b = dict(L=Lb)

            if bridge:
                mf_boltz = mean_field_stats.calc_boltz_binding_cnf_bridge
            else:
                mf_boltz = mean_field_stats.calc_boltz_binding_cnf_loop

            meanfield = (
                mf_boltz(system, tether_type_a, tether_type_b) /
                (mean_field_stats.calc_boltz_exclusion(system,
                                                       tether_type_a) *
                 mean_field_stats.calc_boltz_exclusion(system,
                                                       tether_type_b)))

            self.assertAlmostEqual(explicit, meanfield, 4)

        test_it(20 * nm, 20 * nm, 13 * nm, bridge=True)
        test_it(20 * nm, 20 * nm, 20 * nm, bridge=True)
        test_it(20 * nm, 20 * nm, 23 * nm, bridge=True)
        test_it(20 * nm, 20 * nm, 41 * nm, bridge=True)

        test_it(20 * nm, 11 * nm, 13 * nm, bridge=True)
        test_it(20 * nm, 11 * nm, 20 * nm, bridge=True)
        test_it(20 * nm, 11 * nm, 23 * nm, bridge=True)
        test_it(20 * nm, 11 * nm, 41 * nm, bridge=True)

        test_it(11 * nm, 20 * nm, 13 * nm, loop=True)
        test_it(11 * nm, 20 * nm, 20 * nm, loop=True)
        test_it(11 * nm, 20 * nm, 23 * nm, loop=True)
        test_it(11 * nm, 20 * nm, 41 * nm, loop=True)

        test_it(20 * nm, 20 * nm, 13 * nm, loop=True)
        test_it(20 * nm, 20 * nm, 20 * nm, loop=True)
        test_it(20 * nm, 20 * nm, 23 * nm, loop=True)
        test_it(20 * nm, 20 * nm, 41 * nm, loop=True)

        test_it(20 * nm, 11 * nm, 13 * nm, loop=True)
        test_it(20 * nm, 11 * nm, 20 * nm, loop=True)
        test_it(20 * nm, 11 * nm, 23 * nm, loop=True)
        test_it(20 * nm, 11 * nm, 41 * nm, loop=True)

        test_it(11 * nm, 20 * nm, 13 * nm, loop=True)
        test_it(11 * nm, 20 * nm, 20 * nm, loop=True)
        test_it(11 * nm, 20 * nm, 23 * nm, loop=True)
        test_it(11 * nm, 20 * nm, 41 * nm, loop=True)


# Rods grafted on spheres
class TestRodsGraftedOnSpheres(unittest.TestCase):
    """Unit test for RodsGraftedOnSpheres."""

    def test_segment_point_distance(self):
        v1 = np.array([0.0, 0.0, 0.0])
        v2 = np.array([1.0, 1.0, 1.0])

        self.assertAlmostEqual(
            _segment_point_distance2(v1, v2, np.array([0.5, 0.5, 0.5])),
            0.0)
        self.assertAlmostEqual(
            _segment_point_distance2(v1, v2, np.array([1.0, 0.0, 0.0])),
            (1 * sin(acos(1 / sqrt(3)))) ** 2)
        self.assertAlmostEqual(
            _segment_point_distance2(v1, v2, np.array([0.0, 1.0, 0.0])),
            (1 * sin(acos(1 / sqrt(3)))) ** 2)
        self.assertAlmostEqual(
            _segment_point_distance2(v1, v2, np.array([0.0, 0.0, 1.0])),
            (1 * sin(acos(1 / sqrt(3)))) ** 2)

        r = np.array([-1.0, -3.0, -7.0])
        self.assertAlmostEqual(np.sum((r - v1) ** 2),
                               _segment_point_distance2(v1, v2, r))

        r = np.array([+1.0, +3.0, +7.0])
        self.assertAlmostEqual(np.sum((r - v2) ** 2),
                               _segment_point_distance2(v1, v2, r))

    def test_exclusion(self):
        class FakeSystem(object):
            def __init__(self):
                self.periodic = False
                self.sphere_info = dict(A=dict(centre=[0.0, 0.0, 0.0],
                                               radius=5.0),
                                        B=dict(centre=[20.0, 0.0, 0.0],
                                               radius=8.0),
                                        )

        system = FakeSystem()

        tether_i = dict(sphere='A', L=1.0, pos=[5.0, 0.0, 0.0])
        tether_j = dict(sphere='B', L=1.0, pos=[-8.0, 0.0, 0.0])

        stats = RodsGraftedOnSpheres()

        self.assertAlmostEqual(
            1.0,
            stats.calc_boltz_exclusion(system, tether_i))
        self.assertAlmostEqual(
            1.0,
            stats.calc_boltz_exclusion(system, tether_j))

        system.sphere_info['B']['centre'] = [5.0 + 8.0, 0.0, 0.0]
        self.assertAlmostEqual(
            0.0,
            stats.calc_boltz_exclusion(system, tether_i))
        self.assertAlmostEqual(
            0.0,
            stats.calc_boltz_exclusion(system, tether_j))

    def test_binding(self):
        # For now, just test a restricted geometry where the
        # number "good_fraction" in calc_boltz_binding_cnf is always 1.0
        # (i.e., all points on the circle of possible binding points
        # are allowed)
        class FakeSystem(object):
            def __init__(self):
                self.sphere_info = dict(A=dict(centre=[0.0, 0.0, 0.0],
                                               radius=5.0),
                                        B=dict(centre=[20.0, 0.0, 0.0],
                                               radius=8.0),
                                        )

        system = FakeSystem()

        tether_i = dict(sphere='A', L=1.0, pos=[5.0, 0.0, 0.0])
        tether_j = dict(sphere='B', L=1.0, pos=[-8.0, 0.0, 0.0])

        stats = RodsGraftedOnSpheres()

        self.assertEqual(
            0.0,
            stats.calc_boltz_binding_cnf(system, tether_i, tether_j))

        for d in (1.7, 0.4):

            system.sphere_info['B']['centre'] = [5.0 + 8.0 + d, 0.0, 0.0]

            Li = tether_i['L']
            Lj = tether_j['L']
            d_i_to_bind_plane = (d ** 2 + Li ** 2 - Lj ** 2) / (2 * d)
            gamma_i = acos((d ** 2 + Li ** 2 - Lj ** 2) / (2 * d * Li))
            gamma_j = acos((d ** 2 + Li ** 2 - Lj ** 2) / (2 * d * Li))
            Rcirc = Li * sin(gamma_i)  # == Lj*sin(gamma_j)

            omega_i = (2 * pi * Li ** 2 *
                       stats.calc_boltz_exclusion(system, tether_i))
            omega_j = (2 * pi * Lj ** 2 *
                       stats.calc_boltz_exclusion(system, tether_j))
            omega_ij = 2 * pi * Rcirc / sin(gamma_i + gamma_j)

            self.assertAlmostEqual(
                omega_ij,
                stats._omega_bound(system, tether_i, tether_j))

            self.assertAlmostEqual(
                omega_i,
                stats._omega_unbound(system, tether_i))
            self.assertAlmostEqual(
                omega_j,
                stats._omega_unbound(system, tether_j))

    def test_immutable(self):
        # Unit test for squashing bug where _omega_bound changed the
        # position of the tethers (!!!)

        class FakeSystem(object):
            def __init__(self):
                self.sphere_info = dict(A=dict(centre=[2.0, 0.0, 0.0],
                                               radius=5.0),
                                        )
        system = FakeSystem()

        tether_i = dict(sphere='A', L=1.0, pos=np.array((5.0, 0.0, 0.0)))
        tether_j = dict(sphere='A', L=1.0, pos=np.array((0.0, 5.0, 0.0)))

        stats = RodsGraftedOnSpheres()
        stats._omega_bound(system, tether_i, tether_j)
        self.assertEqual(tether_i['pos'][0], 5.0)
        self.assertEqual(tether_j['pos'][1], 5.0)
