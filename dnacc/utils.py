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
==================================================
Utility classes and functions (:mod:`dnacc.utils`)
==================================================

.. rubric:: Classes

.. autosummary::
   :toctree:

   SymDict
   default_zero_dict
   SyntheticList
   csr_matrix_from_dict

.. rubric:: Functions

.. autosummary::
   :toctree:

   csr_matrix_items
   is_csr_matrix_symmetric
   random_point_sphere
   pbc_delta
"""

__all__ = ['SymDict',
           'default_zero_dict',
           'SyntheticList',
           'csr_matrix_from_dict',
           'csr_matrix_items',
           'is_csr_matrix_symmetric',
           'random_point_sphere',
           'pbc_delta',
           ]

import numpy as np
from math import pi, sin, cos, acos, sqrt

# This neat class comes directly from a stackoverflow.com answer by
# Justin Peel (http://stackoverflow.com/questions/4368423 ...
# ... /python-symmetric-dictionary-where-dab-dba)
#
# I've added the __delitem__ and __contains__ methods
#class SymDict(dict):
#    def __getitem__(self, key):
#        return dict.__getitem__(self,
#                                key if key[0] < key[1]
#                                else (key[1],key[0]))
#
#    def __setitem__(self, key, value):
#        dict.__setitem__(self,
#                         key if key[0] < key[1]
#                         else (key[1],key[0]), value)
#
#    def __delitem__(self, key):
#        return dict.__delitem__(self,
#                                key if key[0] < key[1]
#                                else (key[1],key[0]))
#
#    def __contains__(self, key):
#        return dict.__contains__(self,
#                                 key if key[0] < key[1]
#                                 else (key[1],key[0]))


# Actually, it's much faster to just add the items twice
# and forgo the additional indirection when getting items
class SymDict(dict):
    """
    A dictionary for 2-tuples whose order is irrelevant.

    Use exactly as a dict.
    """
    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        if key[0] != key[1]:
            dict.__setitem__(self, (key[1], key[0]), value)

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        if key[0] != key[1]:
            dict.__delitem__(self, (key[1], key[0]))


class default_zero_dict(dict):
    """A dictionary where missing keys have the value 0.0"""
    def __missing__(self, key):
        return 0.0


class SyntheticList(object):
    """A list that calls methods to get, set and delete items."""

    def __init__(self, fgetlen=None, fgetitem=None,
                 fsetitem=None, fdelitem=None):
        self.fgetlen = fgetlen
        self.fgetitem = fgetitem
        self.fsetitem = fsetitem
        self.fdelitem = fdelitem

    def __len__(self):
        if self.fgetlen is None:
            raise TypeError("Can't calculate list length")
        else:
            return self.fgetlen()

    def __getitem__(self, i):
        if self.fgetitem is None:
            raise AttributeError("Can't get item's value")
        else:
            return self.fgetitem(i)

    def __setitem__(self, i, val):
        if self.fsetitem is None:
            raise AttributeError("Can't set item's value")
        else:
            self.fsetitem(i, val)

    def __delitem__(self, i):
        if self.fdelitem is None:
            raise AttributeError("Can't delete item")
        else:
            self.fdelitem(i)


class csr_matrix_from_dict(object):
    """
    Bare-bones implementation of a CSR-format sparse matrix.

    Notes
    -----
    I wrote this to avoid SciPy's horrendous performance when creating
    a CSR matrix from a dictionary of keys.  See documentation for
    csr_matrix for a discussion of the shape, indptr, indices and data
    attributes.

    In DNACC, we only use sparse matrices for keeping track of values, as
    all the arithmetic is done in the C extension module.

    For iterating through matrix items, use csr_matrix_items().

    Parameters
    ----------
    shape : 2-tuple
        a tuple (M, N) specifying the shape of the matrix
    d : dictionary
        a dictionary mapping (row, col) tuples to values

    Attributes
    ----------
    shape
    indptr
    indices
    data

    Examples
    --------
    >>> mtx = csr_matrix_from_dict((3, 3), {(1, 2): 3, (0, 0): 4})
    >>> mtx[1, 2]
    >>> for (i, j), v in csr_matrix_items(mtx): print( i, j, v )
    >>> for (i, j), v in csr_matrix_items(mtx, row=1): print( i, j, v )
    """

    def __init__(self, shape, d):
        # CSR matrices are laid out so that column indices for row i are
        # stored in indices[indptr[i]:indptr[i+1]] and their corresponding
        # values are stored in data[indptr[i]:indptr[i+1]].
        M, N = shape
        nnz_in_row = np.zeros(M, dtype='int32')
        for (i, j) in d:
            nnz_in_row[i] += 1
        nnz = len(d)

        indptr = np.zeros(M + 1, dtype='int32')
        for i in range(M):
            indptr[i + 1] = indptr[i] + nnz_in_row[i]

        indices = np.zeros(nnz, dtype='int32')
        data = np.zeros(nnz)

        next_indptr = np.array(indptr, copy=True)
        for (i, j), v in d.items():
            assert 0 <= i < M
            assert 0 <= j < N
            ind = next_indptr[i]
            indices[ind] = j
            data[ind] = v
            next_indptr[i] += 1

        self.shape = shape
        self.indptr = indptr
        self.indices = indices
        self.data = data

    def __getitem__(self, key):
        """Retrive item (i,j) from the matrix.  VERY SLOW."""
        i, j = key
        M, N = self.shape
        if not (0 <= i < M) or not (0 <= j < N):
            raise KeyError('Item (%d,%d) out of range (matrix is %dx%d)' %
                           (i, j, M, N))

        for ind in range(self.indptr[i], self.indptr[i + 1]):
            if self.indices[ind] == j:
                return self.data[ind]
        else:
            return 0.0


def csr_matrix_items(mtx, row=None):
    """A generator to efficiently iterate through CSR matrix items.

    If the matrix were stored as a dictionary ``d``, this method is
    analogous to ``d.iteritems()``.

    Parameters
    ----------
    mtx : sparse matrix in CSR format
        The matrix over whose items to iterate
    row : integer or None
        If not None, only iterate over the items in this row

    Examples
    --------
    >>> for (i, j), v in csr_matrix_items(mtx): print(i, j, v)
    >>> for (i, j), v in csr_matrix_items(mtx, row=3): print(i, j, v)
    """

    # CSR matrices are laid out so that column indices for row i are
    # stored in indices[indptr[i]:indptr[i+1]] and their corresponding
    # values are stored in data[indptr[i]:indptr[i+1]].
    indices = mtx.indices
    indptr = mtx.indptr
    data = mtx.data

    for i in range(mtx.shape[0]) if row is None else (row,):
        for ind in range(indptr[i], indptr[i + 1]):
            yield (i, indices[ind]), data[ind]


def is_csr_matrix_symmetric(mtx):
    """Test whether or not the matrix mtx is symmetric.

    Notes
    -----
    VERY SLOW, use only in during development.
    """
    for (i, j), Wij in csr_matrix_items(mtx):
        if i != j and mtx[j, i] != Wij:
            return False
    else:
        return True


def random_point_sphere(radius, centre=(0.0, 0.0, 0.0)):
    """Generate a random point (x,y,z) on a sphere

    Parameters
    ----------
    radius : float
        Sphere radius
    centre : 3-tuple of floats
        Sphere centre

    Returns
    -------
    pt : 3-tuple of float
        A random point chosen uniformly over the sphere surface.
    """
    phi = 2 * pi * np.random.rand()
    costheta = 2 * np.random.rand() - 1.0
    theta = acos(costheta)
    x = centre[0] + radius * sin(theta) * cos(phi)
    y = centre[1] + radius * sin(theta) * sin(phi)
    z = centre[2] + radius * cos(theta)
    return (x, y, z)


def pbc_delta(x1, x2, Lx):
    """Find x1 - x2 in a periodic 1D interval of size Lx."""
    dx = (x1 - x2) % Lx
    if dx > 0.5 * Lx:
        dx -= Lx
    return dx
