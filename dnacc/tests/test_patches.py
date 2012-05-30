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
from dnacc.plates import *
from dnacc.spheres import *
from dnacc.patches import *
from dnacc.units import nm
import numpy as np
import matplotlib.pyplot as plt


## Unit testing
## ------------
def test_plates():
    plates = Plates(200 * nm, 200 * nm)

    plates.beta_DeltaG0["alpha", "alpha'"] = -10

    plates.set_tether_prototype(plate='lower', L=20 * nm,
                                sticky_end="alpha")

    for x, y in [(10, 20), (15, 20), (15, 15)]:
        plates.add_tether(pos=(x * nm, y * nm))

    plates.set_tether_prototype(plate='upper', L=20 * nm,
                                sticky_end="alpha'")

    for x, y in [(12, 21), (17, 21), (14, 16)]:
        plates.add_tether(pos=(x * nm, y * nm))

    plates.separation = 10 * nm

    plates.update()
