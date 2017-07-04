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
from dnacc.plates_meanfield import *
import dnacc.units
from nose.tools import *


## Unit testing
## ------------
class TestPlatesMeanField(unittest.TestCase):
    def test_basics(self):
        nm = dnacc.units.nm
        plates = PlatesMeanField()
        plates.separation = 10 * nm

        plates.set_tether_type_prototype(L=20 * nm)
        plates.add_tether_type(plate='lower', sticky_end='alpha',
                               sigma=1 / (20 * nm) ** 2)
        plates.add_tether_type(plate='upper', sticky_end='alphap',
                               sigma=1 / (20 * nm) ** 2)
        plates.beta_DeltaG0['alpha', 'alphap'] = -10

        plates.update()

        self.assertAlmostEqual(plates.free_energy_density,
                               -0.00109129204, 5)

        self.assertAlmostEqual(plates.binding_free_energy_density,
                               -0.00455702794, 5)

        self.assertAlmostEqual(plates.rep_free_energy_density,
                               0.003465735902, 5)

        exp_sigma_free = [0.00070127133076655418,
                          0.00070127133076655418]
        for exp, got in zip(exp_sigma_free, plates.sigma_free):
            self.assertAlmostEqual(got, exp)

        exp_sigma_bound = [[0, 0.00179873],
                           [0.00179873, 0]]

        for i in range(0, len(plates.tether_types)):
            for j in range(0, len(plates.tether_types)):
                self.assertAlmostEqual(plates.sigma_bound[i, j],
                                       exp_sigma_bound[i][j])

    def test_find_tether_types(self):
        nm = dnacc.units.nm

        plates = PlatesMeanField()
        plates.separation = 10 * nm

        plates.set_tether_type_prototype(sigma=1 / (20 * nm) ** 2)
        i1 = plates.add_tether_type(
            plate='lower', L=20 * nm, sticky_end='alpha')
        i2 = plates.add_tether_type(
            plate='lower', L=10 * nm, sticky_end='alpha')
        i3 = plates.add_tether_type(
            plate='upper', L=20 * nm, sticky_end='alpha')
        i4 = plates.add_tether_type(
            plate='upper', L=20 * nm, sticky_end='beta')

        self.assertRaises(ValueError, plates.find_tether_types)

        self.assertEqual(plates.find_tether_types(sticky_end='gamma'),
                          set())
        self.assertEqual(plates.find_tether_types(sticky_end='alpha'),
                          set((i1, i2, i3)))
        self.assertEqual(plates.find_tether_types(sticky_end='beta'),
                          set((i4,)))
        self.assertEqual(plates.find_tether_types(plate='lower'),
                          set((i1, i2)))
        self.assertEqual(plates.find_tether_types(plate='upper'),
                          set((i3, i4)))
        self.assertEqual(plates.find_tether_types(plate='upper',
                                                   sticky_end='alpha'),
                          set((i3,)))
        self.assertEqual(plates.find_tether_types(plate='lower',
                                                   sticky_end='alpha'),
                          set((i1, i2)))
        self.assertEqual(plates.find_tether_types(L=10 * nm),
                          set((i2,)))
        self.assertEqual(plates.find_tether_types(L=10 * nm,
                                                   sticky_end='beta'),
                          set())
