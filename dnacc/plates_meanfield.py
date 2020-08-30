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

__all__ = ['PlatesMeanField']

from .generic import Generic
from .utils import SymDict, SyntheticList, csr_matrix_from_dict
from math import exp, log
from . import units
from .tether_statistics import RodsGraftedOnPlatesMeanField
import numpy as np


class PlatesMeanField(object):
    """Two parallel DNA-coated plates, treated at a mean-field level.

    Parameters
    ----------
    tether_statistics : object/class
       an object or class that specifies the polymer statistics for the
       grafted tethers.  See the :mod:`.tether_statistics` module for
       details.  If set to None, use a :class:`RodsGraftedOnPlatesMeanField`
       object.

    Attributes
    ----------
    tether_types
        tether_types[a] is a dictionary of the properties of a-type tethers,
        used by the tether_statistics object to calculate free energies of
        binding and exclusion.

    beta_DeltaG0
        beta_DeltaG0[alpha, beta] is the solution hybridisation energies of
        sticky ends alpha and beta.  The labels correspond to the
        ``sticky_end`` properties of individual tether types (see
        :func:`add_tether_type`).

    separation
        Perpendicular distance between plates

    p_free
    sigma_free
    p_bound
    sigma_bound
    sigma_bonds
    binding_free_energy_density
    rep_free_energy_density
    free_energy_density


    See also
    --------
    Plates : explicit-tether treatment of parallel plates
    """

    def __init__(self, tether_statistics=None):
        self.tether_statistics = (tether_statistics if tether_statistics
                                  else RodsGraftedOnPlatesMeanField())
        self._tether_type_prototype = {}
        self.tether_types = []
        self.beta_DeltaG0 = SymDict()

        # Plate-specific info
        self.separation = 0.0

        # Raw free energies in reference state
        self._ref_binding_free_energy_density = 0.0
        self._ref_rep_free_energy_density = 0.0

        # Internal Generic object that drives the calculation
        self._dnacc = None

        # p_free and p_bound here means something slightly different that
        # in a Generic, see p_bound property docstring
        self._p_free = SyntheticList(
            fgetitem=lambda i: (self._dnacc.p_free[i] /
                                self._dnacc.weights[i]))
        self._p_bound = SyntheticList(
            fgetitem=lambda pair: (self._dnacc.p_bound[pair] /
                                   self._dnacc.weights[pair[0]]))

    def set_tether_type_prototype(self, **kwargs):
        """Set properties common to all tether types added after this call.

        The properties common to all new tethers are passed in as keyword
        arguments.

        See :func:`add_tether_type` for details
        """
        self._tether_type_prototype = dict(kwargs)

    def add_tether_type(self, **kwargs):
        """Add a new tether type to this system.

        The keyword arguments specify the tether type's properties, on top
        of those of the prototype tether type (see
        :func:`set_tether_type_prototype`).  Tether type properties are used
        by the tether_statistics object to calculate free energies of
        binding and exclusion.

        Each tether type should have at least three attributes:

          'plate'
            identifies plate on which this tether type is grafted.  Two
            tethers are on the same/opposite plates if their ``plate``
            properties are equal/different.

          'sigma'
            area density of grafting points

          'sticky_end'
            identity of sticky end of tethers of this type, used to look up
            solution hybridisation energy in beta_DeltaG0.

        Different tether_statistics objects may require additional
        properties to be specified.

        Returns
        -------
        idx : int
            Index of newly-added tether type

        Examples
        --------
        >>> plates = PlatesMeanField()
        >>> plates.set_tether_type_prototype(plate='lower', L=20 * nm)
        >>> plates.add_tether_type(sigma=1 / (20 * nm) ** 2,
                                   sticky_end='alpha')
        >>> plates.add_tether_type(sigma=2 / (20 * nm) ** 2,
                                   sticky_end='beta')
        >>> plates.add_tether_type(plate='upper',
                                   sigma=2 / (20 * nm) ** 2,
                                   sticky_end='gamma')
        >>> plates.beta_DeltaG0['alpha', 'beta'] = -2
        """

        new_type_info = self._tether_type_prototype.copy()
        new_type_info.update(kwargs)

        # Check that the attributes 'plate' and 'sigma' are defined
        # and have sensible values
        if not 'plate' in new_type_info:
            raise ValueError("Tether type must have a plate attribute")

        if not 'sigma' in new_type_info:
            raise ValueError("Tether type must have a sigma attribute")

        # Other checks specific to particular tether statistics (e.g.,
        # check rod length is defined for rigid rods)
        try:
            self.tether_statistics.check_tether_type(self, new_type_info)
        except AttributeError:
            pass

        self.tether_types.append(new_type_info)
        return len(self.tether_types) - 1

    def find_tether_types(self, **kwargs):
        """Find tether types with the given properties.

        The properties to match are specified as keyword arguments.

        Returns
        -------
        match_set : set of int
            A set with the indices of the matching tether types.

        Examples
        --------
        Find all tether types on the 'upper' plate

        >>> upper_set = plates.find_tether_types(plate='upper')

        Find all tether types on the 'upper' plate with sticky end 'alpha'

        >>> upper_alpha_set = plates.find_tether_types(plate='upper',
                                                       sticky_end='alpha')
        """

        if not kwargs:
            raise ValueError("No properties specified in find_tether_types")

        result = set()
        for i, t in enumerate(self.tether_types):

            if (all(prop in t for prop in kwargs) and
                all(t[prop] == val for prop, val in kwargs.items())):

                result.add(i)

        return result

    def at(self, h):
        """Place the system at separation h.

        Parameters
        ----------
        h : float
            New separation between plates.

        Returns
        -------
        self : :class:`PlatesMeanField` object
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

    def set_reference_now(self):
        """Define the zero of free energy as the current free energy."""

        self._ref_binding_free_energy_density = 0.0
        self._ref_rep_free_energy_density = 0.0

        self._ref_binding_free_energy_density = \
            self.binding_free_energy_density
        self._ref_rep_free_energy_density = self.rep_free_energy_density

        assert abs(self.binding_free_energy_density < 1e-8)
        assert abs(self.rep_free_energy_density < 1e-8)

    def update(self):
        """
        Update internals following changes in tether types and beta_DeltaG0.
        """

        N = len(self.tether_types)
        boltz_bind = SymDict()

        stats = self.tether_statistics
        loop = stats.calc_boltz_binding_cnf_loop
        bridge = stats.calc_boltz_binding_cnf_bridge

        # Run checks on system and tethers
        try:
            stats.check_system(self)
        except AttributeError:
            pass

        try:
            for t in self.tether_types:
                stats.check_tether_type(self, t)
        except AttributeError:
            pass

        # Set up weights as grafting densities
        weights = np.array(list(t['sigma'] for t in self.tether_types),
                           dtype='float64')

        # Set up repulsion
        boltz_rep = [stats.calc_boltz_exclusion(self, t)
                     for t in self.tether_types]
        self._rep_free_energy_density = (
            -sum(t['sigma'] * log(x)
                 for (t, x) in zip(self.tether_types, boltz_rep)))

        # Set up binding
        for i, t_i in enumerate(self.tether_types):
            for j, t_j in enumerate(self.tether_types):

                if i <= j:
                    continue

                binding_pair = (t_i['sticky_end'], t_j['sticky_end'])

                if binding_pair in self.beta_DeltaG0:

                    if t_i['plate'] == t_j['plate']:
                        boltz_cnf = loop(self, t_i, t_j)
                    else:
                        boltz_cnf = bridge(self, t_i, t_j)
                    boltz_cnf /= boltz_rep[i] * boltz_rep[j]

                    boltz_soln = exp(-self.beta_DeltaG0[binding_pair])
                    boltz_bind[i, j] = boltz_soln * boltz_cnf

        # Calculate p_free, p_bound and avg_num_bound
        self._dnacc = Generic(csr_matrix_from_dict((N, N), boltz_bind),
                              weights)

    def _dnacc_check(self):
        if self._dnacc is None:
            raise RuntimeError("You must call update() whenever you "
                               "change any tether properties or binding "
                               "strengths")

    @property
    def p_free(self):
        """p_free[a] = probability that an a-type tether is unbound."""
        self._dnacc_check()
        return self._p_free

    @property
    def sigma_free(self):
        """sigma_free[a] = area density of unbound a-type tethers."""
        self._dnacc_check()
        return self._dnacc.p_free

    @property
    def p_bound(self):
        """
        p_bound[a, b] = probability that one a-type tether is bound """\
        """to any b-type tether.
        """
        self._dnacc_check()
        return self._p_bound

    @property
    def sigma_bound(self):
        """
        sigma_bound[a, b] = area density of bonds between a-type """\
        """and b-type tethers.
        """
        self._dnacc_check()
        return self._dnacc.p_bound

    @property
    def sigma_bonds(self):
        """Area density of bonds."""
        self._dnacc_check()
        return self._dnacc.avg_num_bonds

    @property
    def binding_free_energy_density(self):
        """Binding free energy per unit area of system.

        This excludes repulsion due to volume exclusion."""
        self._dnacc_check()
        return (self._dnacc.binding_free_energy -
                self._ref_binding_free_energy_density)

    @property
    def rep_free_energy_density(self):
        """Free energy of repulsion per unit area due to volume exclusion."""
        self._dnacc_check()
        return (self._rep_free_energy_density -
                self._ref_rep_free_energy_density)

    @property
    def free_energy_density(self):
        """Free energy per unit area of system."""
        return (self.rep_free_energy_density +
                self.binding_free_energy_density)
