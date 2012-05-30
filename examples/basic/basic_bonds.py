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

plates.beta_DeltaG0['alpha', 'alphap'] = -10  # in kT

h_arr = np.linspace(1 * nm, 40 * nm, 40)
num_bonds_arr = [plates.at(h).avg_num_bonds for h in h_arr]

print("# h (nm)     Bonds Formed / Max possible")
for (h, n) in zip(h_arr, num_bonds_arr):
    print (h / nm), (n / (num_tethers / 2))
