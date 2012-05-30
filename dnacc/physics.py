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
==========================================================================
Physical constants in :mod:`correct units <units>` (:mod:`dnacc.physics``)
==========================================================================

Fundamental constants
+++++++++++++++++++++

.. autodata:: c
.. autodata:: mu_0
.. autodata:: eps_0
.. autodata:: N_A
.. autodata:: kB
.. autodata:: R

Atomic and Nuclear Physics
++++++++++++++++++++++++++

.. autodata:: e
.. autodata:: eV
.. autodata:: amu
.. autodata:: m_e
.. autodata:: a_0
.. autodata:: D
"""

from . import units

_GSL_CONST_MKSA_SPEED_OF_LIGHT = 2.99792458e8  # m / s
_GSL_CONST_MKSA_VACUUM_PERMITTIVITY = 8.854187817e-12  # A^2 s^4 / kg m^3
_GSL_CONST_MKSA_VACUUM_PERMEABILITY = 1.25663706144e-6  # kg m / A^2 s^2
_GSL_CONST_NUM_AVOGADRO = 6.02214199e23  # 1 / mol
_GSL_CONST_MKSA_BOLTZMANN = 1.3806504e-23  # kg m^2 / K s^2
_GSL_CONST_MKSA_MOLAR_GAS = 8.314472e0  # kg m^2 / K mol s^2
_GSL_CONST_MKSA_ELECTRON_CHARGE = 1.602176487e-19  # A s
_GSL_CONST_MKSA_ELECTRON_VOLT = 1.602176487e-19  # kg m^2 / s^2
_GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS = 1.660538782e-27  # kg
_GSL_CONST_MKSA_MASS_ELECTRON = 9.10938188e-31  # kg
_GSL_CONST_MKSA_BOHR_RADIUS = 5.291772083e-11  # m
_GSL_CONST_MKSA_DEBYE = 3.33564095198e-30  # A s^2 / m^2

# Fundamental constants
# =====================

#: Speed of light in vacuum
c = _GSL_CONST_MKSA_SPEED_OF_LIGHT * units.m / units.s

#: Permeability of free space, :math:`\mu_0`
mu_0 = _GSL_CONST_MKSA_VACUUM_PERMEABILITY * units.N / units.Ampere ** 2

#: Permittivity of free space, :math:`\epsilon_0`
eps_0 = _GSL_CONST_MKSA_VACUUM_PERMITTIVITY * units.F / units.m

#: Avogadro's number
N_A = _GSL_CONST_NUM_AVOGADRO

#: Boltzmann's constant
kB = _GSL_CONST_MKSA_BOLTZMANN * units.J / units.K

#: Gas constant (numerically identical to kB in base units!)
R = _GSL_CONST_MKSA_MOLAR_GAS * units.J / (units.K * units.mol)


# Atomic and Nuclear Physics
# ==========================

#: Electron charge
e = _GSL_CONST_MKSA_ELECTRON_CHARGE * units.C

#: Electron volt
eV = _GSL_CONST_MKSA_ELECTRON_VOLT * units.J

#: Atomic mass unit
amu = _GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS * units.kg

#: Mass of electron
m_e = _GSL_CONST_MKSA_MASS_ELECTRON * units.kg

#: Bohr radius
a_0 = _GSL_CONST_MKSA_BOHR_RADIUS * units.m

#: Debye
D = _GSL_CONST_MKSA_DEBYE * units.C * units.m
