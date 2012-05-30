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

__all__ = ['Spheres']

from .system import System
from .tether_statistics import RodsGraftedOnSpheres


class Spheres(System):
    """A collection of DNA-coated spheres, with explicit tethers.

    :class:`Spheres` inherits all the attributes and methods of
    :class:`System`.

    Parameters
    ----------
    tether_statistics : object/class
       an object or class that specifies the polymer statistics for the
       grafted tethers.  See the :mod:`.tether_statistics` module for
       details.  If set to None, use a :class:`RodsGraftedOnSpheres` object.

    Attributes
    ----------
    sphere_info
        A dictionary of spheres and their properties (see
        :func:`add_sphere`)

    Methods
    -------
    add_sphere
    sphere_centre_at

    Examples
    --------
    Create two spheres, one with radius 50 nm, the other, 100 nm, coated
    with rigid-rod-like DNA tethers. Note that tether positions are
    three-dimensional vectors from the sphere centres

    >>> from dnacc.units import nm
    >>> spheres = Spheres()
    >>> spheres.add_sphere('A', (0.0, 0.0, 0.0), 50.0 * nm)
    >>> spheres.add_sphere('B', (200.0 * nm, 0.0, 0.0), 100.0 * nm)
    >>> spheres.add_tether(sphere='A', pos=(50 * nm,0,0), L=10 * nm)
    >>> spheres.add_tether(sphere='B', pos=(-100 * nm,0,0), L=10 * nm)
    """

    def __init__(self, tether_statistics=None):
        super(Spheres, self).__init__(tether_statistics if tether_statistics
                                      else RodsGraftedOnSpheres())

        self.sphere_info = {}

    def add_sphere(self, name, centre, radius):
        """Add a sphere to the system.

        Parameters
        ----------
        name : object
            An identifier (e.g. a string, a number) used to refer to the
            new sphere

        centre : 3-tuple of floats
            Coordinates of sphere centre

        radius : float
            Sphere radius
        """
        self.sphere_info[name] = dict(centre=centre, radius=radius)

    def sphere_centre_at(self, name, centre, update=True):
        """Move a sphere to a new position.

        Parameters
        ----------
        name : object
            Identifies sphere to move (see :func:`add_sphere`)

        centre : 3-tuple of floats
            Coordinates of new sphere centre

        update : bool
            Whether or not to call :func:`update` after moving the sphere

        Returns
        -------
        self : :class:`Spheres` object
            self object, for ready access to an attribute (see example).

        Examples
        --------
        >>> spheres.sphere_centre_at('red_ball', (0,0,0)).free_energy

        Notes
        -----
        If you use this method to move many spheres, set update to False
        for all of the calls, then call :func:`update` at the end.
        """

        if name not in self.sphere_info:
            raise ValueError('No sphere "%s"' % str(name))

        self.sphere_info[name]['centre'] = centre

        if update:
            self.update()

        return self
