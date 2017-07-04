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
=======================================================================
Consistent units definitions for molecular systems (:mod:`dnacc.units`)
=======================================================================

.. currentmodule:: dnacc.units

The internal base units are GROMACS units, which are useful,
consistent, and internally yield numbers that are close to 1 in
most atomistic setups:

- Length = nm
- Mass = amu
- Time = ps (together yield natural energy unit of kJ/mol)
- Current = e / ps
- Temperature = K

Electrostatic formulas have all the usual SI prefactors.  We don't
bother with the units relating to light (cd and related units) since
they are mostly irrelevant for condensed matter simulations.

Initially, we express the SI base units in terms of GROMACS units,
then derive all other units from them.  For internal convenience,
we actually use g instead of kg for defining masses.

Almost all unit definitions also set up the related prefixed units
from femto-X to tera-X.

Whenever you read in a unitful quantity from the user, multiply
it by the relevant units.  For example,

   >>> myLength = fields[5] * units.nm

Whenever you output unitful quantities, divide by the units you want
to use.  For example,

   >>> print("Total energy = %.2f kcal/mol" %
             (E / (units.kcal / units.mol)))

.. note::

  The unit "Celsius" is not defined explicitly to not confuse it with
  "Coulombs".

Base units
++++++++++

All of these can be prefixed with the usual SI prefixes, e.g. ``nm`` for
nanometer.

.. autodata:: m
.. autodata:: g
.. autodata:: s
.. autodata:: Ampere
.. autodata:: K
.. autodata:: mol

By default, the unit for Amperes is :data:`.Ampere`, not ``A`` as usual.
To change this behaviour, call :func:`.add_amperes_unit`.

.. autofunction:: add_amperes_unit

Derived units
+++++++++++++

Most of these can also be prefixed with the usual SI prefixes, e.g. ``pN``
for pico-Newton.

.. autodata:: Hz
.. autodata:: rad
.. autodata:: sr
.. autodata:: N
.. autodata:: Pa
.. autodata:: J
.. autodata:: W
.. autodata:: C
.. autodata:: V
.. autodata:: F
.. autodata:: Ohm
.. autodata:: S
.. autodata:: Wb
.. autodata:: T
.. autodata:: H

Non-SI units
++++++++++++

.. autodata:: min
.. autodata:: h
.. autodata:: d
.. autodata:: degree
.. autodata:: arcmin
.. autodata:: arcsec
.. autodata:: ha
.. autodata:: L
.. autodata:: t

Physics units
+++++++++++++

Some units which are more properly considered physical constants are
defined in the :mod:`.physics` module.

.. autodata:: G
.. autodata:: bar
.. autodata:: atm
.. autodata:: Torr
.. autodata:: mmHg
.. autodata:: P

Chemistry units
+++++++++++++++

Units that pop up regularly in chemistry.

.. autodata:: AA
.. autodata:: cal
.. autodata:: kcal
.. autodata:: M
.. autodata:: cc

"""

from math import pi

_GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS = 1.660538782e-27  # kg
_GSL_CONST_MKSA_ELECTRON_CHARGE = 1.602176487e-19  # A s
_GSL_CONST_NUM_AVOGADRO = 6.02214199e23  # 1 / mol
_GSL_CONST_MKSA_GAUSS = 1e-4  # kg / A s^2
_GSL_CONST_MKSA_BAR = 1e5  # kg / m s^2
_GSL_CONST_MKSA_STD_ATMOSPHERE = 1.01325e5  # kg / m s^2
_GSL_CONST_MKSA_TORR = 1.33322368421e2  # kg / m s^2
_GSL_CONST_MKSA_METER_OF_MERCURY = 1.33322368421e5  # kg / m s^2

#: Meter
m = 1e9
#: Gram
g = 1e-3 / _GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS
#: Second
s = 1e12
#: Ampere
Ampere = (1. / _GSL_CONST_MKSA_ELECTRON_CHARGE) / s
#: Kelvin
K = 1.0
#: mole
mol = _GSL_CONST_NUM_AVOGADRO


def add_amperes_unit():
    """Define ``A`` as the unit Ampere.

    By default, this definition is disabled to avoid accidental mixups with
    the unit for Angstroms (:data:`.AA`).  This mixup can result in some
    thoroughly puzzling bugs."""

    global A
    A = Ampere


_SI_prefixes = {
    'f': 1e-15,
    'p': 1e-12,
    'n': 1e-9,
    'u': 1e-6,
    'm': 1e-3,
    'c': 1e-2,
    'd': 1e-1,
    'da': 1e+1,
    'h': 1e+2,
    'k': 1e+3,
    'M': 1e+6,
    'G': 1e+9,
    'T': 1e+12
}


def _add_prefixes(name, realname=None):
    """Add standard SI prefixes to a base unit."""
    globs = globals()
    val = globs[realname or name]
    for prefix, factor in list(_SI_prefixes.items()):
        globs[prefix + name] = factor * val

_add_prefixes('m')
_add_prefixes('g')
_add_prefixes('s')
_add_prefixes('A', 'Ampere')
_add_prefixes('K')
_add_prefixes('mol')

# Named derived units, after
# http://en.wikipedia.org/wiki/SI_derived_unit

#: Hertz
Hz = 1 / s
#: Radian
rad = m / m
#: Steradian
sr = m ** 2 / m ** 2
#: Newton
N = kg * m / s ** 2
#: Pascal
Pa = N / m ** 2
#: Joule
J = N * m
#: Watt
W = J / s
#: Coulomb
C = Ampere * s
#: Volt
V = J / C
#: Farad
F = C / V
#: Ohm
Ohm = V / Ampere
#: Siemens
S = 1 / Ohm
#: Weber
Wb = J / Ampere
#: Tesla
T = Wb / m ** 2
#: Henry
H = Wb / Ampere
# degree C is dangerous, confused with Coulombs
#degC = K                # Celsius
# Some of the other units are irrelevant for what I do...

_add_prefixes('Hz')
_add_prefixes('N')
_add_prefixes('Pa')
_add_prefixes('J')
_add_prefixes('W')
_add_prefixes('C')
_add_prefixes('V')
_add_prefixes('F')
_add_prefixes('Ohm')

# Now some official non-SI units, after
# http://en.wikipedia.org/wiki/Non-SI_units_accepted_for_use_with_SI

#: minute
min = 60 * s
#: hour
h = 60 * min
#: day
d = 24 * h

#: degree of arc
degree = (pi / 180.0)
#: arc-minute
arcmin = degree / 60.0
#: arc-second
arcsec = arcmin / 60.0

#: hectare
ha = 10000 * m ** 2

#: litre
L = dm ** 3

#: tonne
t = 1e3 * kg

# Some physics-based units are defined in physics.py, some here

#: Gauss
G = _GSL_CONST_MKSA_GAUSS * T

#: bar (prefixed versions available, e.g. mbar = millibar)
bar = _GSL_CONST_MKSA_BAR * Pa
_add_prefixes('bar')
#: atmosphere
atm = _GSL_CONST_MKSA_STD_ATMOSPHERE * Pa
#: Torr
Torr = _GSL_CONST_MKSA_TORR * Pa
#: mmHg
mmHg = 1e-3 * _GSL_CONST_MKSA_METER_OF_MERCURY * Pa

# Viscosity
#: Poise (prefixed versions available, e.g. cP = centipoise)
P = 1 * g / (cm * s)
_add_prefixes('P')

# Some useful units in chemistry
# ==============================

#: Angstrom.  See also :func:`.add_amperes_unit`.
AA = 1e-10 * m

cal = 4.184 * J
"""Thermochemical calorie

.. note::

  The calorie here is the thermochemical calorie (4.184 J), not the
  International Steam Table calorie (4.1868 J) used in GSL.
"""

#: kilocalorie
kcal = 1e3 * cal

#: molar (prefixed versions available, e.g. mM = millimolar)
M = mol / L
_add_prefixes('M')

#: cubic centimeter
cc = cm ** 3
