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

"""
Derjaguin approximation sphere-sphere potentials from plate-plate potentials.
"""

__all__ = ['calc_spheres_potential']

from math import pi
import numpy as np
from scipy.integrate import cumtrapz


def calc_spheres_potential(h_arr, V_plate_arr, R1, R2=0.0):
    """Approximate sphere-sphere potential from plate-plate potential.

    Use the Derjaguin approximation to calculate the interaction potential
    between two spheres as a function of separation, starting from the
    interaction potential per unit area of two uniformly coated plates.

    Parameters
    ----------
    h_arr : list of floats
      a list of separations at which the plate-plate interaction strength
      per unit area has been calculated.

    V_plate_arr : list of floats
      V_plate_arr[i] is the plate-plate interaction strength per unit area
      at separation h_arr[i]

    R1, R2 : float
      radii of spheres.  If R2 is 0.0, take R2 = R1  If R2 >> R1, the
      results is effectively a sphere-plate potential.

    Returns
    -------
    V_sphere_arr : list of floats
      V_sphere_arr[i] = the sphere-sphere interaction when the distance of
      closest approach between the spheres is hArr[i].

    Notes
    -----
    The Derjaguin approximation takes the following form:

    .. math::

       V_{sphere}(h) = \\frac{2 \\pi}{R_1^{-1} + R_2^{-1}}
                         \\int_h^\\infty dh'\\,V_{plate}(h')

    It results from assuming that each patch on one sphere interacts with
    the closest point on the other sphere, and then integrating this
    interaction per unit area over the surface of the sphere.  The
    approximation is effective when both R1 and R2 are far larger than the
    range of the plate-plate potential.
    """

    if R2 == 0.0:
        R2 = R1

    prefactor = 2 * pi / (1 / R1 + 1 / R2)

    # Integral done in the reverse direction
    cumint = np.concatenate(([0.0], cumtrapz(V_plate_arr, h_arr)))

    # Reverse and set starting point to zero
    cumint = cumint[-1] - cumint

    return prefactor * cumint
