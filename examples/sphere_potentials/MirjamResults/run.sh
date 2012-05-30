#!/bin/bash

#****************************************************************************************************
#****************************************************************************************************
#*		Copyright: Mirjam E. Leunissen, FOM Institute AMOLF, The Netherlands - 11 October 2010		*
#*		No restrictions apply to the usage of this program provided that users cite:				*
#*		M.E. Leunissen and D. Frenkel, J. Chem. Phys. 134, 084702 (2011), doi: 10.1063/1.3557794	*
#****************************************************************************************************
#****************************************************************************************************

#****************************************************************************************************	
#		This program generates the DNA-mediated pair interaction potentials for planar surfaces		*
#		and spherical particles, functionalized with stiff DNA constructs of Ldna = 20 nm long at	*
#		the user-defined average strand spacing or surface coverage (S), the DNA binding			*
#		strength (G) and the particle radius (R).													*
#		The derivation of the analytical expressions and the working range of these pair			*
#		interactions can be found in Section IV and the rest of:									*
#		M.E. Leunissen and D. Frenkel, J. Chem. Phys. 134, 084702 (2011), doi: 10.1063/1.3557794	*
#****************************************************************************************************

#****************************************************************************************************
#		Set the following parameters below:															*
#		M: mode 1 = calculate plate-plate interaction; 2 = calculate sphere-sphere interaction		*
#		Hmin: minimum surface-surface separation (> 0, in units of Ldna)							*
#		Hmax: maximum surface-surface separation (<= 2, in units of Ldna)							*
#		dH: step size of the calculated data points between Hmin and Hmax							*
#																									*
#		S: average strand spacing (in units of Ldna)												*
#		G: the hybridization free energy of the DNA sticky ends when free in solution (in kT)		*
#		R: particle radius (in units of Ldna)														*
#****************************************************************************************************

M=(1 2)
Hmin=(0.05)
Hmax=(2.0)
dH=(0.05)

S=(0.25 0.75)
G=(-2.0 -3.0 -4.0 -5.0 -6.0 -7.0)
R=(6.7 25.0)

#****************************************************************************************************
#****************************************************************************************************

mkdir ./results

for m in ${M[*]}; do
for s in ${S[*]}; do
for g in ${G[*]}; do
for r in ${R[*]}; do

cat > input <<EOF

M  Hmin  Hmax   dH  S   G   R
$m $Hmin $Hmax $dH  $s  $g  $r    
 
EOF

./DNAPairPotentials

rm input

if [ $m -eq 1 ] ; then
mv potential ./results/plates-S$s-G$g.dat
else
mv potential ./results/spheres-R$r-S$s-G$g.dat
fi

done
done
done
done

echo
echo Copyright: Mirjam E. Leunissen, FOM Institute AMOLF, The Netherlands
echo All users cite: M.E. Leunissen and D. Frenkel, J. Chem. Phys. 134, 084702, 2011, doi: 10.1063/1.3557794
echo
