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

__all__ = ['Plates']

from .system import System
from .tether_statistics import RodsGraftedOnPlates


class Plates(System):
    """Two parallel DNA-coated plates, with explicit tethers.

    :class:`Plates` inherits all the attributes and methods of
    :class:`System`.

    Parameters
    ----------
    Lx, Ly : float
        plate dimensions in the x and y directions (the z direction is
        perpendicular to the plates).

    tether_statistics : object/class
       an object or class that specifies the polymer statistics for the
       grafted tethers.  See the :mod:`.tether_statistics` module for
       details.  If set to None, use a :class:`RodsGraftedOnPlates` object.

    periodic : bool
        whether or not to impose periodic boundary conditions along x and y.

    Attributes
    ----------
    dims : 2-tuple of floats
        Plate dimensions in x and y directions.

    periodic : bool
        Whether or not the system is periodic in x and y.

    separation : float
        Perpendicular distance between plates

    Methods
    -------
    at

    See also
    --------
    PlatesMeanField : mean field treatment of parallel plates

    Examples
    --------
    Create two periodic 200x200 nm plates coated with rigid-rod-like DNA
    tethers.

    >>> plates = Plates( 200 * units.nm, 200 * units.nm )
    """

    def __init__(self, Lx, Ly, tether_statistics=None, periodic=True):
        super(Plates, self).__init__(tether_statistics if tether_statistics
                                     else RodsGraftedOnPlates())
        self.dims = (Lx, Ly)
        self.periodic = periodic
        self.separation = 0

    def at(self, h):
        """Place the system at separation h.

        Parameters
        ----------
        h : float
            New separation between plates.

        Returns
        -------
        self : :class:`Plates` object
            self object, for ready access to an attribute (see example).

        Examples
        --------
        >>> from dnacc.units import nm
        >>> print(plates.at(10 * nm).binding_free_energy)
        """

        if h != self.separation:
            self.separation = h
            self.update()

        return self
