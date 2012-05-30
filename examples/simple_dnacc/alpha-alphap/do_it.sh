#!/bin/bash

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

../../../simple_dnacc CONTROL_BONDS-h-20nm-explicit.txt
../../../simple_dnacc CONTROL_BONDS-h-20nm-meanfield.txt

../../../simple_dnacc CONTROL_POTENTIAL_dg-4kT-explicit.txt
../../../simple_dnacc CONTROL_POTENTIAL_dg-14kT-explicit.txt
../../../simple_dnacc CONTROL_POTENTIAL_dg-24kT-explicit.txt

../../../simple_dnacc CONTROL_POTENTIAL_dg-4kT-meanfield.txt
../../../simple_dnacc CONTROL_POTENTIAL_dg-14kT-meanfield.txt
../../../simple_dnacc CONTROL_POTENTIAL_dg-24kT-meanfield.txt

gnuplot plot.gp
