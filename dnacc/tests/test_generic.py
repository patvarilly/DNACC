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
from dnacc.utils import SymDict, csr_matrix_from_dict
import numpy as np
from math import *
from dnacc.generic import *


## Unit testing
## ------------

# DNACC logic
class TestGeneric(unittest.TestCase):
    """Unit test for DNACC logic."""

    def test_simple_dnacc(self):
        cases = [

            dict(beta_DeltaG0=-10, weights=None,
                 binding_free_energy=-9.01,
                 p_bound_0_1=0.99327740159843414),

            dict(beta_DeltaG0=-10, weights=[0.8, 0.6, 3.7],
                 binding_free_energy=-5.54,
                 p_bound_0_1=0.8 * 1.2497165073687115 * 0.6),

            dict(beta_DeltaG0=-1, weights=None,
                 binding_free_energy=-1.05,
                 p_bound_0_1=0.55013119441689895),

            dict(beta_DeltaG0=-1, weights=[0.8, 0.6, 3.7],
                 binding_free_energy=-0.59,
                 p_bound_0_1=0.8 * 0.69788894233575782 * 0.6)]

        for case in cases:
            N = 3
            boltz_bind = SymDict()
            boltz_bind[0, 1] = exp(-case['beta_DeltaG0'])
            boltz_bind = csr_matrix_from_dict((N, N), boltz_bind)

            dnacc = Generic(boltz_bind, case['weights'])

            self.assertAlmostEqual(dnacc.binding_free_energy,
                                   case['binding_free_energy'], 2)
            self.assertAlmostEqual(dnacc.p_bound[0, 1],
                                   case['p_bound_0_1'])
            self.assertAlmostEqual(dnacc.p_bound[0, 1],
                                   dnacc.p_bound[1, 0])
            self.assertEqual(dnacc.p_bound[0, 0], 0.0)
            self.assertEqual(dnacc.p_bound[1, 2], 0.0)
            self.assertTrue(dnacc.avg_num_bonds >= 0)

    def test_count_bonds(self):
        N = 3
        boltz_bind = SymDict()
        boltz_bind[0, 1] = boltz_bind[0, 2] = boltz_bind[1, 2] = 100.0
        boltz_bind = csr_matrix_from_dict((N, N), boltz_bind)

        dnacc = Generic(boltz_bind)

        self.assertAlmostEqual(dnacc.count_bonds(set((0,)), set((1,))),
                               0.46587246970021012)
        self.assertAlmostEqual(dnacc.count_bonds(set(range(N)), set(range(N))),
                               dnacc.avg_num_bonds, delta=1e-5)
