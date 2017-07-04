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

__all__ = ['System']

from .generic import Generic
from .utils import SymDict, csr_matrix_from_dict
from math import exp, log
from . import units
from . import tether_statistics
import numpy as np


class System(object):
    """A system of DNA-coated colloids with explicit tethers.

    Parameters
    ----------
    tether_statistics : object/class
       an object or class that specifies the polymer statistics for the
       grafted tethers.  See the :mod:`.tether_statistics` module for
       details.  A commonly used tether_statistics is
       :class:`.RodsGraftedOnPlates`.

    Attributes
    ----------
    tethers
        tethers[i] is a dictionary of the properties of tether i, used by
        the tether_statistics object to calculate free energies of binding
        and exclusion.

    beta_DeltaG0
        beta_DeltaG0[alpha, beta] is the solution hybridisation energies of
        sticky ends alpha and beta.  The labels correspond to the
        ``sticky_end`` properties of individual tethers (see
        :func:`add_tether`).

    num_tethers
    p_free
    p_bound
    avg_num_bonds
    binding_free_energy
    rep_free_energy
    free_energy

    Notes
    -----
    A useful specialization is :class:`.Plates`.  A mean-field
    approximation of it is :class:`.PlatesMeanField`.

    After creating a :class:`System` object and populating it with tethers
    (:func:`add_tether`), you **must** call :func:`update` to calculate the
    values of all the system properties.  If you subsequently change the
    system, you must call :func:`update` again.  A fast update is possible
    if the only change was the values of :attr:`beta_DeltaG0`.

    """

    def __init__(self, tether_statistics):
        self.tether_statistics = tether_statistics
        self.tethers = []
        self.beta_DeltaG0 = SymDict()

        # Raw free energies in reference state
        self._ref_binding_free_energy = 0.0
        self._ref_rep_free_energy = 0.0

        # A dictionary with default tether properties for new tethers
        # built with :func:`add_tether`.
        self._tether_prototype = {}

        # Configurational Boltzmann factors (useful when we only
        # need to update binding energies)
        self._boltz_binding_cnf = SymDict()

        # Internal Generic object that drives the calculation
        self._dnacc = None

    def set_tether_prototype(self, **kwargs):
        """Set properties common to all tethers added after this call.

        The properties common to all new tethers are passed in as keyword
        arguments.

        See :func:`add_tether` for details
        """
        self._tether_prototype = dict(kwargs)

    def add_tether(self, **kwargs):
        """Add a new tether to this system.

        The keyword arguments specify the tether's properties, on top of
        those of the prototype tether (see :func:`set_tether_prototype`).
        Tether properties are used by the tether_statistics object to
        calculate free energies of binding and exclusion.

        Each tether should have at least one attribute:

          'sticky_end'
            identity of sticky end of tether, used to look up solution
            hybridisation energy in beta_DeltaG0.

        Different tether_statistics objects may require additional
        properties to be specified.

        Returns
        -------
        idx : int
            Index of newly-added tether

        Examples
        --------

        >>> system.set_tether_prototype(plate='lower', L=20 * nm,
                                        sticky_end="alpha")
        >>> system.add_tether(pos=(10*nm, 20*nm))
        >>> system.add_tether(pos=(15*nm, 20*nm))
        >>> system.add_tether(pos=(15*nm, 15*nm))
        >>> system.set_tether_prototype(plate='upper', L=20 * nm,
                                        sticky_end="alpha'")
        >>> system.add_tether(pos=(12*nm, 21*nm))
        >>> system.add_tether(pos=(17*nm, 21*nm))
        >>> system.add_tether(pos=(14*nm, 16*nm))
        """
        new_tether_info = self._tether_prototype.copy()
        new_tether_info.update(kwargs)

        # Check that sticky_end attribute is defined
        if not 'sticky_end' in new_tether_info:
            raise ValueError('Must specify "sticky_end" attribute of tether')

        # Other checks specific to particular tether statistics (e.g.,
        # check rod length is defined for rigid rods)
        try:
            self.tether_statistics.check_tether(self, new_tether_info)
        except AttributeError:
            pass

        self.tethers.append(new_tether_info)
        return len(self.tethers) - 1

    def find_tethers(self, **kwargs):
        """Find tethers with the given properties.

        The properties to match are specified as keyword arguments.

        Returns
        -------
        match_set : set of int
            A set with the indices of the matching tethers.

        Examples
        --------
        Find all tethers on the 'upper' plate

        >>> upper_set = plates.find_tethers(plate='upper')

        Find all tethers on the 'upper' plate with sticky end 'alpha'

        >>> upper_alpha_set = plates.find_tethers(plate='upper',
                                                  sticky_end='alpha')
        """

        if not kwargs:
            raise ValueError("No properties specified in find_tethers")

        result = set()
        for i, t in enumerate(self.tethers):

            if (all(prop in t for prop in kwargs) and
                all(t[prop] == val for prop, val in kwargs.items())):

                result.add(i)

        return result

    def update(self, DeltaG0_only=False):
        """
        Update internals following changes in tethers and beta_DeltaG0.

        Parameters
        ----------
        DeltaG0_only : bool
           Set to True if only ``beta_DeltaG0`` has changed since the last
           update.  This avoids recalculating configurational binding
           entropies, whose cost is quite significant.

        Notes
        -----
        If previously, ``beta_DeltaG0`` was infinite and now it is not, that
        is considered a geometry change.  In other words, if the set of
        tethers that can bind can change, you **must** set DeltaG0_only to
        False.
        """

        N = len(self.tethers)
        boltz_bind = SymDict()

        stats = self.tether_statistics
        calc_boltz_binding_cnf = stats.calc_boltz_binding_cnf

        # Run checks on system and tethers
        try:
            stats.check_system(self)
        except AttributeError:
            pass

        try:
            for t in self.tethers:
                stats.check_tether(self, t)
        except AttributeError:
            pass

        # Set up repulsion
        boltz_rep = [stats.calc_boltz_exclusion(self, t)
                     for t in self.tethers]
        self._rep_free_energy = -sum(log(x) for x in boltz_rep)

        # Pre-exponentiate beta_DeltaG0
        boltz_soln = SymDict()
        for (i, j), beta_DeltaG0 in self.beta_DeltaG0.items():
            if i <= j:
                boltz_soln[i, j] = exp(-beta_DeltaG0)
        ends = [t['sticky_end'] for t in self.tethers]

        # Do we need to recalulate configurational Boltzmann factors?
        if not DeltaG0_only:

            self._boltz_binding_cnf.clear()

            # Enumerate interaction partners of each strand
            if hasattr(stats, 'calc_interaction_lists'):
                interaction_list = stats.calc_interaction_lists(self)
            else:
                N = len(self.tethers)
                interaction_list = [range(i + 1, N) for i in range(N)]

            # Set up binding
            for i, nns_i in enumerate(interaction_list):
                info_i = self.tethers[i]

                for j in nns_i:
                    if i > j:
                        continue
                    info_j = self.tethers[j]

                    sticky_end_pair = (ends[i], ends[j])
                    if sticky_end_pair not in boltz_soln:
                        continue

                    boltz_cnf = (
                        calc_boltz_binding_cnf(self, info_i, info_j) /
                        (boltz_rep[i] * boltz_rep[j]))
                    if boltz_cnf != 0:
                        self._boltz_binding_cnf[i, j] = boltz_cnf
                        boltz_bind[i, j] = (boltz_soln[sticky_end_pair] *
                                            boltz_cnf)

        else:

            # Simply run through old _boltz_binding_cnf factors and
            # multiply by new exp(-beta_DeltaG0) factors
            for (i, j), boltzCnf in self._boltz_binding_cnf.items():
                boltz_bind[i, j] = boltz_soln[ends[i], ends[j]] * boltzCnf

        # Calculate p_free, p_bound and avg_num_bound
        self._dnacc = Generic(csr_matrix_from_dict((N, N), boltz_bind))

    def set_reference_now(self):
        """Define the zero of free energy as the current free energy."""

        self._ref_binding_free_energy = 0.0
        self._ref_rep_free_energy = 0.0

        self._ref_binding_free_energy = self.binding_free_energy
        self._ref_rep_free_energy = self.rep_free_energy

        assert abs(self.binding_free_energy < 1e-8)
        assert abs(self.rep_free_energy < 1e-8)

    def _dnacc_check(self):
        if self._dnacc is None:
            raise RuntimeError("You must call update() whenever "
                               "you change any tether properties or "
                               "binding strengths")

    @property
    def num_tethers(self):
        """Number of tethers in the system."""
        return len(self.tethers)

    @property
    def p_free(self):
        """p_free[i] = probability that tether i is unbound."""
        self._dnacc_check()
        return self._dnacc.p_free

    @property
    def p_bound(self):
        """p_bound[i,j] = probability that tethers i and j are bound."""
        self._dnacc_check()
        return self._dnacc.p_bound

    @property
    def avg_num_bonds(self):
        """Average number of bonds."""
        self._dnacc_check()
        return self._dnacc.avg_num_bonds

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
        self._dnacc_check()
        return self._dnacc.count_bonds(i_set, j_set)

    @property
    def binding_free_energy(self):
        """Binding free energy of system.

        This excludes repulsion due to volume exclusion."""
        self._dnacc_check()
        return (self._dnacc.binding_free_energy -
                self._ref_binding_free_energy)

    @property
    def rep_free_energy(self):
        """Free energy of repulsion due to volume exclusion."""
        self._dnacc_check()
        return self._rep_free_energy - self._ref_rep_free_energy

    @property
    def free_energy(self):
        """Free energy of system."""
        return self.binding_free_energy + self.rep_free_energy
