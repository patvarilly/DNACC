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
from dnacc.physics import *


## Unit testing
## ------------
class TestPhysics(unittest.TestCase):
    def test_basics(self):
        self.assertAlmostEqual(e / units.ps, 1.0)

    def test_thermal(self):
        kT = kB * 300 * units.K
        self.assertAlmostEqual(kT / (units.kJ / units.mol), 2.49, 2)
        self.assertAlmostEqual(kT / (units.kcal / units.mol), 0.60, 2)

    def test_electrostatic(self):
        from math import pi
        self.assertAlmostEqual(D / (e * units.AA), 0.21, 2)
        self.assertAlmostEqual(e ** 2 / units.AA / (4 * pi * eps_0) /
                               (units.kcal / units.mol), 332.06372)

    def test_water_properties(self):
        rho_w = 0.997 * units.g / units.cc
        m_w = 18.015268 * amu
        self.assertAlmostEqual((rho_w / m_w / units.M), 55, 0)

    def test_quantum(self):
        from math import pi
        E_0 = 0.5 * (1 / (4 * pi * eps_0)) * (e ** 2 / a_0)
        self.assertAlmostEqual(E_0 / eV, 13.6, 1)
