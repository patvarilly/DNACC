The directory hdist_H.NN contains the Boltzmann binding factors and derived
data for two spheres separated by H Kuhn lengths (5 nm).  NN refers to a
particular realisation of the random tether positions on the two spheres (there
should be 15 of these at each separation).

In each directory, tethers.dat gives the 3D absolute coordinates (in units
of Kuhn lengths) of the grafting points of a and b tethers on the spheres.
The left sphere is centered at (-R,0,0) and the right sphere is centered at
(R+h,0,0).

The file hyb1.dat lists all non-zero values of exp(-beta Delta G_ij) as
follows: there are N_A blocks of data, corresponding to the N_A a-type
strands.  In each block of data, there are n_{aB} lines, where n_{aB} is the
number of b-type strands that this a strand can bind to.  Note that the
strand indices are independent, i.e., strands a go from 1 to N_A, and
strands b go from 1 to N_B.

hyb2.dat contains identical information, but from the point of view of the B
strands.

rep1.dat is a list of the partition function of the confined vs unconfined a
tethers, while rep2.dat is the analogous list for the b tethers.
