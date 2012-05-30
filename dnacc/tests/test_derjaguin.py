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
from dnacc.derjaguin import *
import numpy as np
from math import pi


## Unit testing
## ------------
class TestDerjaguin(unittest.TestCase):
    def test_simple(self):
        # Try a square well potential
        depth = 3
        hArr = np.linspace(0, 10.0, 100)
        VArr = -depth * np.ones(100)

        # Expect a linear potential with separation
        R = 1000.0
        VSphereArrExp = pi * R * -depth * (hArr[-1] - hArr)
        VSphereArrGot = calc_spheres_potential(hArr, VArr, R)

        for exp, got in zip(VSphereArrExp, VSphereArrGot):
            self.assertAlmostEqual(exp, got)
