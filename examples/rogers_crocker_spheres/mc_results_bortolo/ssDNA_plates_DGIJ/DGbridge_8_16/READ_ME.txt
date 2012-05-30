

configurational entropy cost as a function of the planes distance (h) and
the angle (theta_i) made between the tether-tether direction and the normal
plane vector


Notes on the choosen unit: in Q_ij we do not have yet divided by rho_0 (5 nm)^3  (=75.2768d0). This has been done in the program that compute n_hyb vs DGa


=================================================================

FILE STRUCTURE (notice the first two lines)

Number of tethered computed
dtheta
 ..................
theta_i   ; Q_ij/Q_i Q_j (or exp(-\beta DGconf))
 ..................

theta_i goes from 0 to ACos[h/16] - dtheta

==================================================================

