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
A generic DNA-coated colloid (DNACC) system.

At an abstract level, a DNACC system consists of a set of N tethers that can
bind pairwise with binding free energy :math:`\\Delta G_{ij}`.  This binding
energy may be infinite if i and j physically cannot bind.  Otherwise, it
consists of the solution hybridisation free energy of the sticky end
sequences of i and j along with the configuration free energy cost of
binding.
"""

__all__ = ['Generic', '_do_calculate_py', '_do_add_up_py']

# I've modelled the extension code on Michael Shirts' and John Chodera's
# "pymbar" Python implementation of their MBAR algorithm
# (see https://simtk.org/home/pymbar)

import numpy as np
import scipy.integrate
from math import log, exp
from .utils import (csr_matrix_items, csr_matrix_from_dict,
                    default_zero_dict, is_csr_matrix_symmetric)

# Try to import C extension module to speed the inner loops here
try:
    from . import _generic
    _use_C_extension = True
except ImportError:
    _use_C_extension = False
    print("Could not import extension code to speed up inner loops\n"
          "-- using much slower pure Python version instead.")


class Generic(object):
    """A generic DNA-coated colloid (DNACC) system.

    Everywhere, N is the number of tethers in the DNACC system.  This object
    is immutable, i.e., once the parameters are specified at creation, they
    cannot be changed later.

    Parameters
    ----------
    boltz_bind : 2D sparse matrix in CSR format
      an NxN symmetric matrix, with boltz_bind[i,j] = :math:`K_{ij} =
      \\exp(-\\beta \\Delta G_{ij})`, i.e. the ratio of partition functions
      of the tethers i and j when bound vs. unbound.  The easiest way to
      create this matrix is by first building a dictionary of entries, then
      passing it to :func:`.csr_matrix_from_dict`.

    weights : None, or list of N floats
      If not None, a weight to assign to each tether's probability of being
      free (see Notes for details).  This is used by
      :class:`.PlatesMeanField` to take into account the grafting densities
      of each type of tether.  The calculated :attr:`p_free` and
      :attr:`p_bound` are also weighted.

    Other Parameters
    ----------------
    self_consistent_max_delta : float
      Threshold change in any :math:`p_i` to terminate fixed-point
      iteration.
    self_consistent_max_steps : integer
      Maximum number of fixed-point iterations.
    thermo_int_epsabs : float
      Absolute maximum error, in kT, allowed in thermodynamic integration
      for computing binding free energy.
    thermo_int_epsrel : float
      Absolute relative error allowed in thermodynamic integration
      for computing binding free energy.

    Attributes
    ----------
    avg_num_bonds
      average number of bonds formed
    p_free
      p_free[i] = probability that tether i is unbound.
    p_bound
      p_bound[i,j] = probability that tethers i and j are bound.

    binding_free_energy

    Notes
    -----
    If weights are not specified, solve the following self-consistent
    equations for the probability of tether i being unbound, :math:`p_i`:

    .. math:: p_i = 1 / (1 + \sum_j K_{ij} p_j),

    whereby the probability of i and j being bound is:

    .. math:: p_{ij} = p_i K_{ij} p_j.

    If weights are specified, then solve the following self-consistent
    equations for the weighted probability of tether i being unbound,
    :math:`s_i = w_i p_i`:

    .. math:: s_i = w_i / (1 + \sum_j K_{ij} s_j)

    whereby the weighted probability of i and j being bound is:

    .. math:: s_{ij} = s_i K_{ij} s_j.

    The accuracy of the fixed point iteration can be controlled using the
    attributes ``self_consistent_max_delta`` and
    ``self_consistent_max_steps``.

    The binding energy is computed using the following thermodynamic
    integration:

    .. math::
       F_{att} = - \\int_0^{10\\,k_B T - \\min(\\Delta G_{ij})}
                           d\\lambda\, \\sum_{i < j} p_{ij}(\\lambda)

    where :math:`p_{ij}(\\lambda)` is the probability that tethers i and j
    are bound after shifting all binding energies :math:`\\Delta G_{ij}` up
    by an amount :math:`\\lambda`.  The integration is performed to absolute
    accuracy ``thermo_int_epsabs`` and relative accuracy
    ``thermo_int_epsrel``.
    """

    def __init__(self, boltz_bind, weights=None,
                 self_consistent_max_delta=1e-7,
                 self_consistent_max_steps=10001,
                 thermo_int_epsabs=1e-2,
                 thermo_int_epsrel=1e-3):

        # Description of the system
        M, N = boltz_bind.shape
        assert M == N
        #assert is_csr_matrix_symmetric( boltz_bind )

        self._boltz_bind = boltz_bind
        self._weights = None if weights is None else np.array(weights)

        # Variables controlling solution accuracy
        self._self_consistent_max_delta = self_consistent_max_delta
        self._self_consistent_max_steps = self_consistent_max_steps
        self._thermo_int_epsabs = thermo_int_epsabs
        self._thermo_int_epsrel = thermo_int_epsrel

        # Space for results
        self.avg_num_bonds = -1
        self.p_free = np.zeros(N)
        self.p_bound = default_zero_dict()

        # Run self-consistent calculation
        self._calculate()

        # Postpone thermodynamic integration
        self._thermo_int_ready = False

    def _calculate(self, boltz_prefactor=1.0):
        """Recalculate self-consistent values of p_free and avg_num_bonds.
        """

        # Main calculation loop
        mtx = self._boltz_bind

        if _use_C_extension:

            delta, step = _generic.do_calculate(
                self.p_free, mtx.indptr, mtx.indices, mtx.data,
                boltz_prefactor,
                self._self_consistent_max_delta,
                self._self_consistent_max_steps,
                self._weights)

            self.avg_num_bonds = (
                boltz_prefactor *
                _generic.do_add_up(self.p_free, mtx.indptr,
                                   mtx.indices, mtx.data))

        else:

            # See replacement Python code below
            delta, step = _do_calculate_py(
                self.p_free, mtx,
                boltz_prefactor,
                self._self_consistent_max_delta,
                self._self_consistent_max_steps,
                self._weights)

            self.avg_num_bonds = (
                boltz_prefactor *
                _do_add_up_py(self.p_free, mtx))

        # Finish off by calculating p_bound and weighted counterparts
        # (skip this part during thermodynamic integration)
        if boltz_prefactor == 1.0:
            # To every non-zero element in _boltz_bind, there's a non-zero
            # element in p_bound.
            self.p_bound.clear()
            p = self.p_free
            for (i, j), Kij in csr_matrix_items(mtx):
                self.p_bound[i, j] = p[i] * Kij * p[j]

    def _thermodynamic_integrand(self, DeltaDeltaG):
        """Return sum_(i<=j) p_ij when
        betaDeltaG[i,j] -> betaDeltaG[i,j] + DeltaDeltaG"""

        self._calculate(boltz_prefactor=exp(-DeltaDeltaG))
        return -self.avg_num_bonds

    def _calc_free_energies(self):

        assert not self._thermo_int_ready

        # Calc binding free energy using thermodynamic integration over
        # DeltaDeltaG, where betaDeltaG[i,j] -> betaDeltaG[i,j] +
        # DeltaDeltaG

        # The upper integration limit is technically infinity, but
        # using a limit whereby the lowest DeltaG is +10 kT is good enough
        # If nothing binds, then clearly the free energy of binding is 0.0
        try:
            # The maximum Boltzmann factor corresponds to the binding
            # with lowest (i.e. strongest) DeltaG0
            max_boltz_bind = max(self._boltz_bind.data)
            maxL = -(-log(max_boltz_bind))
            maxL += 10

        except ValueError:
            maxL = 0.0

        if maxL <= 0.0:
            # Nothing binds, easy peasy
            self._binding_free_energy = 0.0
            self._thermo_int_ready = True
            return

        # Save old values of p_free and avg_num_bonds
        old_p_free = np.array(self.p_free, copy=True)
        old_avg_num_bonds = self.avg_num_bonds

        # Calculate integral with SciPy
        self._binding_free_energy = scipy.integrate.quad(
            self._thermodynamic_integrand, 0, maxL,
            epsabs=self._thermo_int_epsabs,
            epsrel=self._thermo_int_epsrel)[0]

        # Place system back in a sane state
        self.p_free = old_p_free
        self.avg_num_bonds = old_avg_num_bonds

        self._thermo_int_ready = True

    @property
    def binding_free_energy(self):
        """Free energy of binding of this system.

        Calculating this attribute is expensive, so its evaluation is
        postponed until the attribute's value is first requested."""

        if not self._thermo_int_ready:
            self._calc_free_energies()
        return self._binding_free_energy

    def count_bonds(self, i_set, j_set):
        """Counts bonds between tethers in i_set and j_set.

        Parameters
        ----------
        i_set : set of int
            One set of tether ids
        j_set : set of int
            Another set of tether ids

        Notes
        -----
        For performance, it's *essential* that i_set and j_set be sets and
        not lists.
        """

        # The spirit of this code is as follows:
        #count = 0.0
        #for i in i_set:
        #    partial_count = 0.0
        #    for j in j_set:
        #        partial_count += self._boltz_bind[i,j] * self.p_free[j]
        #    count += self.p_free[i] * partial_count
        #return count

        # But for performance, it's usually better to iterate
        # on the items of the _boltz_bind matrix, since that's the
        # sparse entity
        mtx = self._boltz_bind
        p = self.p_free

        count = 0.0
        for i in i_set:
            if p[i] != 0.0:
                partial_count = 0.0
                for (gobble, j), Wij in csr_matrix_items(mtx, row=i):
                    if j in j_set:
                        partial_count += Wij * p[j]
                count += p[i] * partial_count

        return count


# Replacement Python methods, in case C extension module cannot be loaded
def _do_calculate_py(p_free, mtx, prefactor,
                    maxDelta, maxSteps, weights):

    if weights is None:
        if np.any(p_free > 1.0) or any(p_free < 0.0):
            p_free[:] = 1.0
    else:
        if np.any(p_free > weights) or any(p_free < 0.0):
            p_free[:] = weights

    delta = 1.0
    step = 0
    while delta > maxDelta and step < maxSteps:
        delta = 0.0
        for i in xrange(len(p_free)):
            oldP = p_free[i]
            #x = mtx[i, :] * p_free
            x = 0.0
            for (ii, j), Wij in csr_matrix_items(mtx, row=i):
                x += Wij * p_free[j]
            w = weights[i] if weights is not None else 1.0
            p_free[i] = w / (1.0 + prefactor * x)
            delta = max(delta, abs(oldP - p_free[i]) / (w or 1.0))
        step += 1

    return delta, step


def _do_add_up_py(p_free, mtx):
    avgN = 0.0
    for i in xrange(len(p_free)):
        if p_free[i] != 0.0:
            cumSum = 0.0
            for (ii, j), Wij in csr_matrix_items(mtx, row=i):
                if i <= j:
                    cumSum += Wij * p_free[j]
            avgN += p_free[i] * cumSum
    return avgN
