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

import dnacc
from dnacc.units import nm
import numpy as np

num_tethers = 200

box_L = 200 * nm
plates = dnacc.Plates(Lx=box_L, Ly=box_L, periodic=True)

for (x, y) in np.random.random_sample((num_tethers / 2, 2)):
    plates.add_tether(plate='lower',
                      sticky_end='alpha',
                      L=20 * nm,
                      pos=(x * box_L, y * box_L))

for (x, y) in np.random.random_sample((num_tethers / 2, 2)):
    plates.add_tether(plate='upper',
                      sticky_end='alphap',
                      L=20 * nm,
                      pos=(x * box_L, y * box_L))

plates.beta_DeltaG0['alpha', 'alphap'] = 0.0
plates.at(20 * nm)

print("# Delta G_0 (kT)     Bonds Formed / Max possible")
for dg in np.linspace(-20.0, 0.0, 21):
    plates.beta_DeltaG0['alpha', 'alphap'] = dg
    plates.update(DeltaG0_only=True)
    print dg, (plates.avg_num_bonds / (num_tethers / 2))
