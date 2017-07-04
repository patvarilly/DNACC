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
from dnacc.system import *
import dnacc.units
from dnacc.tether_statistics import RodsGraftedOnPlates


## Unit testing
## ------------
class TestSystem(unittest.TestCase):
    def test_basics(self):
        nm = dnacc.units.nm
        plates = System(RodsGraftedOnPlates())
        plates.dims = (200 * nm, 200 * nm)
        plates.periodic = True
        plates.separation = 10 * nm

        plates.beta_DeltaG0["alpha", "alpha'"] = -10

        plates.set_tether_prototype(plate='lower', L=20 * nm,
                                    sticky_end="alpha")

        for x, y in [(10, 20), (15, 20), (15, 15)]:
            plates.add_tether(pos=(x * nm, y * nm))

        plates.set_tether_prototype(plate='upper', L=20 * nm,
                                    sticky_end="alpha'")

        for x, y in [(12, 21), (17, 21), (14, 16)]:
            plates.add_tether(pos=(x * nm, y * nm))

        plates.update()

        self.assertAlmostEqual(plates.free_energy, -1.45, 2)

        self.assertAlmostEqual(plates.binding_free_energy, -5.61, 2)

        self.assertAlmostEqual(plates.rep_free_energy, 4.16, 2)

        exp_p_free = [0.29917599, 0.22935016, 0.29560323,
                      0.25106227, 0.31209903, 0.2609686]
        for exp, got in zip(exp_p_free, plates.p_free):
            self.assertAlmostEqual(got, exp, 6)

        exp_p_bound = [[0, 0, 0, 0.4267047, 0.1334154, 0.1407036],
                       [0, 0, 0, 0.2103240, 0.4066406, 0.1536849],
                       [0, 0, 0, 0.1119088, 0.1478448, 0.4446428],
                       [0.4267047, 0.2103240, 0.1119088, 0, 0, 0],
                       [0.1334154, 0.4066406, 0.1478448, 0, 0, 0],
                       [0.1407036, 0.1536849, 0.4446428, 0, 0, 0]]

        for i in range(0, plates.num_tethers):
            for j in range(0, plates.num_tethers):
                self.assertAlmostEqual(plates.p_bound[i, j],
                                       exp_p_bound[i][j], 6)

    def test_find_tethers(self):
        nm = dnacc.units.nm

        plates = System(RodsGraftedOnPlates())
        plates.dims = (200 * nm, 200 * nm)
        plates.periodic = True
        plates.separation = 10 * nm

        plates.set_tether_prototype(pos=(0.0, 0.0))
        i1 = plates.add_tether(plate='lower', L=20 * nm, sticky_end='alpha')
        i2 = plates.add_tether(plate='lower', L=10 * nm, sticky_end='alpha')
        i3 = plates.add_tether(plate='upper', L=20 * nm, sticky_end='alpha')
        i4 = plates.add_tether(plate='upper', L=20 * nm, sticky_end='beta')

        self.assertRaises(ValueError, plates.find_tethers)

        self.assertEqual(plates.find_tethers(sticky_end='gamma'),
                          set())
        self.assertEqual(plates.find_tethers(sticky_end='alpha'),
                          set((i1, i2, i3)))
        self.assertEqual(plates.find_tethers(sticky_end='beta'),
                          set((i4,)))
        self.assertEqual(plates.find_tethers(plate='lower'),
                          set((i1, i2)))
        self.assertEqual(plates.find_tethers(plate='upper'),
                          set((i3, i4)))
        self.assertEqual(plates.find_tethers(plate='upper',
                                              sticky_end='alpha'),
                          set((i3,)))
        self.assertEqual(plates.find_tethers(plate='lower',
                                              sticky_end='alpha'),
                          set((i1, i2)))
        self.assertEqual(plates.find_tethers(L=10 * nm),
                          set((i2,)))
        self.assertEqual(plates.find_tethers(L=10 * nm,
                                              sticky_end='beta'),
                          set())
