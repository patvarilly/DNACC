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

#Mean Field 
../../../simple_dnacc CONTROL_BONDS_h_10nm_meanfield.txt
grep "3 plate2/4 plate2" BONDS_MF >aa
grep "1 plate1/2 plate1" BONDS_MF >bb
paste aa bb > cc
cat cc | awk '{print $1,$2+$7}'>BONDS_mean_h=1.00-weak
rm aa bb cc

grep "6 plate2" BONDS_MF >BONDS_mean_h=1.00-ggp

grep "6 plate2" BONDS_MF >aa
grep "3 plate2/1 plate1" BONDS_MF >bb
paste aa bb > cc
cat cc | awk '{print $1,$2+$7}'>BONDS_mean_h=1.00-strong
rm aa bb cc

../../../simple_dnacc CONTROL_POTENTIAL_dg-5kT_meanfield.txt
../../../simple_dnacc CONTROL_POTENTIAL_dg-15kT_meanfield.txt
../../../simple_dnacc CONTROL_POTENTIAL_dg-25kT_meanfield.txt


#Explicit 
../../../simple_dnacc CONTROL_BONDS_h_10nm_explicit.txt
grep "((3, 'plate2'), (4, 'plate2'))" BONDS > aa
grep "((1, 'plate1'), (2, 'plate1'))" BONDS > bb
paste aa bb > cc
cat cc | awk '{print $1,$2+$8}'>BONDS_explicit_h=1.00-weak
rm aa bb cc

grep "((5, 'plate1'), (6, 'plate2'))" BONDS >BONDS_explicit_h=1.00-ggp

grep "((5, 'plate1'), (6, 'plate2'))" BONDS >aa
grep "((1, 'plate1'), (3, 'plate2'))" BONDS > bb
paste aa bb > cc
cat cc | awk '{print $1,$2+$8}'>BONDS_explicit_h=1.00-strong
rm aa bb cc

../../../simple_dnacc CONTROL_POTENTIAL_dg-5kT_explicit.txt
../../../simple_dnacc CONTROL_POTENTIAL_dg-15kT_explicit.txt
../../../simple_dnacc CONTROL_POTENTIAL_dg-25kT_explicit.txt

gnuplot plot.gp
