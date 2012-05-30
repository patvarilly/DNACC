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

plates = dnacc.PlatesMeanField()

plates.add_tether_type(plate='lower',
                       sticky_end='alpha',
                       L=20 * nm,
                       sigma=1 / (20 * nm) ** 2)

plates.add_tether_type(plate='upper',
                       sticky_end='alphap',
                       L=20 * nm,
                       sigma=1 / (20 * nm) ** 2)

plates.beta_DeltaG0['alpha', 'alphap'] = -8  # in kT

plates.at(41 * nm).set_reference_now()

h_arr = np.linspace(1 * nm, 40 * nm, 40)
V_plate_arr = [plates.at(h).free_energy_density for h in h_arr]
R = 500 * nm
V_sphere_arr = dnacc.calc_spheres_potential(h_arr, V_plate_arr, R)

print("# h (nm)     V (kT)")
for (h, V) in zip(h_arr, V_sphere_arr):
    print (h / nm), V
