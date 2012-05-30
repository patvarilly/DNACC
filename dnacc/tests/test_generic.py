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

# C extension module
class TestExtenstion(unittest.TestCase):
    """Unit test for C extension module that speeds up inner loops."""

    def setUp(self):
        self.p_free = np.arange(5, 10).astype('float64')
        d = SymDict()
        d[1, 2] = 600.0
        d[2, 3] = 6000.0
        self.mtx = csr_matrix_from_dict((5, 5), d)
        self.prefactor = 32.0
        self.maxDelta = 1e-7
        self.maxSteps = 10001
        self.weights = np.array([5., 6., 7., 8., 9.])

    def test_all(self):
        try:
            from dnacc import _generic

            for weights in [self.weights, None]:
                self.weights = weights
                p_free_copy1 = self.p_free.copy()
                p_free_copy2 = self.p_free.copy()

                deltaExt, stepExt = _generic.do_calculate(
                    p_free_copy1,
                    self.mtx.indptr, self.mtx.indices, self.mtx.data,
                    self.prefactor, self.maxDelta, self.maxSteps,
                    self.weights)

                deltaPy, stepPy = _do_calculate_py(
                    p_free_copy2, self.mtx, self.prefactor,
                    self.maxDelta, self.maxSteps, self.weights)

                self.assertEqual((deltaExt, stepExt),
                                 (deltaPy, stepPy))

                for p1, p2 in zip(p_free_copy1, p_free_copy2):
                    self.assertAlmostEqual(p1, p2)

                self.assertEqual(
                    _generic.do_add_up(p_free_copy1,
                                       self.mtx.indptr,
                                       self.mtx.indices,
                                       self.mtx.data),
                    _do_add_up_py(p_free_copy2, self.mtx))

        except ImportError:
            pass


# DNACC logic
class TestGeneric(unittest.TestCase):
    """Unit test for DNACC logic."""

    def test_simple_dnacc(self):
        global _use_C_extension

        for use_ext in [True, False]:

            cases = [

                dict(beta_DeltaG0=-10, weights=None,
                     binding_free_energy=-9.01,
                     p_bound_0_1=0.99329197089893884),

                dict(beta_DeltaG0=-10, weights=[0.8, 0.6, 3.7],
                     binding_free_energy=-5.54,
                     p_bound_0_1=0.8 * 1.2497165073687115 * 0.6),

                dict(beta_DeltaG0=-1, weights=None,
                     binding_free_energy=-1.05,
                     p_bound_0_1=0.55013119441689895),

                dict(beta_DeltaG0=-1, weights=[0.8, 0.6, 3.7],
                     binding_free_energy=-0.59,
                     p_bound_0_1=0.8 * 0.69788894233575782 * 0.6),

                ]

            for case in cases:
                N = 3
                boltz_bind = SymDict()
                boltz_bind[0, 1] = exp(-case['beta_DeltaG0'])
                boltz_bind = csr_matrix_from_dict((N, N), boltz_bind)

                _use_C_extension = use_ext
                dnacc = Generic(boltz_bind, case['weights'])

                self.assertAlmostEqual(dnacc.binding_free_energy,
                                       case['binding_free_energy'], 2)
                self.assertAlmostEqual(dnacc.p_bound[0, 1],
                                       case['p_bound_0_1'])
                self.assertAlmostEqual(dnacc.p_bound[0, 1],
                                       dnacc.p_bound[1, 0])
                self.assertEqual(dnacc.p_bound[0, 0], 0.0)
                self.assertEqual(dnacc.p_bound[1, 2], 0.0)
                self.assert_(dnacc.avg_num_bonds >= 0)

    def test_count_bonds(self):
        N = 3
        boltz_bind = SymDict()
        boltz_bind[0, 1] = boltz_bind[0, 2] = boltz_bind[1, 2] = 100.0
        boltz_bind = csr_matrix_from_dict((N, N), boltz_bind)

        dnacc = Generic(boltz_bind)

        self.assertAlmostEqual(dnacc.count_bonds(set((0,)), set((1,))),
                               0.4658727174875526)
