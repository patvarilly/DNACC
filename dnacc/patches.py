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
Utility functions to add patches of tethers to spheres and planes.
"""

__all__ = ['stamp_tethers',
           'add_circular_patch_to_plate',
           'add_circular_patch_to_sphere']

import numpy as np
from math import pi, sqrt, cos, sin


def stamp_tethers(system, pts, **kwargs):
    """Graft many identical tether into a system at the given points.

    The 'pos' property of the i-th added tether is set to pts[i].

    Parameters
    ----------
    pts : NxK array of floats
        K-dimensional coordinates of N grafting points
    kwargs : keyword arguments
        Common properties of tethers in patch

    Returns
    -------
    tether_ids : list of integers
        List of ids for the tethers that were added.

    Examples
    --------
    >>> plates = Plates()
    >>> pts = np.random.random_sample((100, 2))
    >>> stamp_tethers(plates, pts, L=20 * nm, sticky_end='alpha',
    >>>               plate='lower')
    """

    ids = list()
    for pt in pts:
        tether_info = dict(kwargs)
        tether_info['pos'] = pt
        id_ = system.add_tether(**tether_info)
        ids.append(id_)

    return ids


def add_circular_patch_to_plate(system, centre, R, N, **kwargs):
    """Add a circular patch of tethers with the given properties to a plate.

    Parameters
    ----------
    system : :class:`.System` object
        The system to which the patch will be added
    centre : 2-tuple of floats
        (x, y) coordinates of the patch centre
    R : float
        Patch radius
    N : integer
        Number of tethers in patch
    kwargs : keyword arguments
        Common properties of tethers in patch

    Returns
    -------
    tether_ids : list of integers
        List of ids for the tethers that were added.

    Examples
    --------
    >>> plates = Plates(200 * nm, 200 * nm, periodic=True)
    >>> add_circular_patch_to_plate(plates, centre=(100 * nm, 100 * nm),
    >>>                             R=30 * nm, N=20,
    >>>                             L=20 * nm, sticky_end='alpha',
    >>>                             plate='lower')
    """

    pts_rel_r2 = np.random.random_sample((N, 1))
    pts_theta = 2 * pi * np.random.random_sample((N, 1))
    pts_r = R * np.sqrt(pts_rel_r2)
    pts = np.hstack((pts_r * np.cos(pts_theta),
                     pts_r * np.sin(pts_theta)))

    return stamp_tethers(system, pts, **kwargs)


def add_circular_patch_to_sphere(system, centre, angular_aperture, N,
                                 **kwargs):
    """Add a circular patch of tethers with the given properties to a sphere.

    The patch is centred at a 3D point ``centre`` relative to a sphere
    centre ``O``.  All tethers are grafted within an angle
    ``angular_aperture`` of the vector ``O`` -> ``centre``.

    Parameters
    ----------
    system : :class:`.System` object
        The system to which the patch will be added
    centre : 3-tuple of floats
        (x, y, z) coordinates of the patch centre, *relative to sphere
        centre*
    angular_aperture : float
        Half-angle of cone (with apex at the sphere centre) that defines
        patch edges.
    N : integer
        Number of tethers in patch
    kwargs : keyword arguments
        Common properties of tethers in patch

    Returns
    -------
    tether_ids : list of integers
        List of ids for the tethers that were added.

    Examples
    --------
    >>> spheres = Spheres()
    >>> spheres.add_sphere(name='red_ball', centre=(0,0,0), radius=500 * nm)
    >>> add_circular_patch_to_sphere(spheres, centre=(500 * nm, 0, 0),
    >>>                              angular_aperture=pi / 6.0, N=100,
    >>>                              L=20 * nm, sticky_end='alpha',
    >>>                              sphere='red_ball')
    """

    # First generate N random points around the north pole of a sphere,
    radius = sqrt(np.dot(centre, centre))

    min_costheta = cos(angular_aperture)
    costheta = (min_costheta +
                (1 - min_costheta) * np.random.random_sample((N, 1)))
    theta = np.arccos(costheta)
    sintheta = np.sin(theta)

    phi = (2 * pi) * np.random.random_sample((N, 1))
    sinphi = np.sin(phi)
    cosphi = np.cos(phi)

    pts = np.hstack((radius * sintheta * cosphi,
                     radius * sintheta * sinphi,
                     radius * costheta))

    # Build rotation matrix to align north pole to O->centre
    e_k, e_i, e_j = np.linalg.qr(
        np.array([centre, [1, 0, 0], [0, 1, 0], [0, 0, 1]]).T)[0].T
    R = np.matrix(np.vstack([e_i, e_j, e_k]).T)

    # Rotate all the points
    pts = np.array(R * pts.T).T

    # And stamp them in
    return stamp_tethers(system, pts, **kwargs)
