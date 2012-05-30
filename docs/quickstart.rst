Quickstart
==========

DNACC is a Python module.  To use it, you have two choices:

* Write a Python script that sets up DNACC objects and drives a calculation.
  If you know how to program and are comfortable with Python, this method is
  powerful and very general.  See :ref:`dnacc-via-python`.

* Use a simplified control script, ``simple_dnacc``, to perform a number of
  simple but common calculations.  If you are not comfortable with Python,
  or are performing a relatively standard calculation, this method is for
  you.  See :ref:`dnacc-via-simple-dnacc`.

.. _dnacc-via-python:

Using DNACC by writing Python scripts
+++++++++++++++++++++++++++++++++++++

After installing DNACC, try the following Python script to compute the pair
potential between two spheres of radius 500 nm uniformly coated with
complementary DNA tethers.  The script first calculates a plate-plate
potential, then using the Derjaguin approximation to derive a sphere-sphere
potential.  DNA tethers are treated (by default) as rigid rods of length L =
20 nm with a sticky end.  The hybridisation free energy of the sticky ends
in solution in this case is -8 kT.  Tethers are grafted at a density of 1
every (20 nm)^2 on two parallel plates, 'upper' and 'lower'.  The zero of
free energy is taken to be the free energy of the system when the
plates/spheres are separated by 41 nm.

::

  import dnacc
  from dnacc.units import nm
  import numpy as np
  
  plates = dnacc.PlatesMeanField()
  
  plates.add_tether_type(plate='lower',
                         sticky_end='alpha',
			 L=20 * nm,
                         sigma=1 / (20 * nm) ** 2)

  plates.add_tether_type(plate='upper',
                         sticky_end='alphap',
			 L=20 * nm,
                         sigma=1 / (20 * nm) ** 2)
			 
  plates.beta_DeltaG0['alpha', 'alphap'] = -8  # in kT

  plates.at(41 * nm).set_reference_now()

  h_arr = np.linspace(1 * nm, 40 * nm, 40)
  V_plate_arr = [plates.at(h).free_energy_density for h in h_arr]
  R = 500 * nm
  V_sphere_arr = dnacc.calc_spheres_potential(h_arr, V_plate_arr, R)

  print("# h (nm)     V (kT)")
  for (h, V) in zip(h_arr, V_sphere_arr):
      print (h / nm), V

You can find this example in the file ``examples/basic/basic_spheres.py``.
You can run it as follows:

::

  cd examples/basic
  python basic_spheres.py > spheres.dat

The result will be stored in spheres.dat.  When plotted (e.g., by running
``gnuplot plot_spheres.gp``), you should get a potential that looks like
this:

.. image:: Spheres.*

Here's a slightly more involved example: a calculation where DNA tethers are
represented explicitly (as opposed to in mean field), and what is output is
the average number of bonds formed per unit area vs. plate separation:

::

  import dnacc
  from dnacc.units import nm
  import numpy as np
  
  num_tethers = 200
  
  box_L = 200 * nm
  plates = dnacc.Plates(Lx=box_L, Ly=box_L, periodic=True)
  
  for (x, y) in np.random.random_sample((num_tethers/2, 2)):
        plates.add_tether(plate='lower',
                          sticky_end='alpha',
                          L=20 * nm,
                          pos=(x * box_L, y * box_L))
  
  for (x, y) in np.random.random_sample((num_tethers/2, 2)):
        plates.add_tether(plate='upper',
                          sticky_end='alphap',
                          L=20 * nm,
                          pos=(x * box_L, y * box_L))
  
  plates.beta_DeltaG0['alpha', 'alphap'] = -10  # in kT
  
  h_arr = np.linspace(1 * nm, 40 * nm, 40)
  num_bonds_arr = [plates.at(h).avg_num_bonds for h in h_arr]
  
  print("# h (nm)     Bonds Formed / Max possible")
  for (h, n) in zip(h_arr, num_bonds_arr):
        print (h / nm), (n / (num_tethers/2))

This example can be found in ``examples/basic/basic_bonds.py``.  Run it by
executing the following commands::

  cd examples/basic
  python basic_bonds.py > bonds.dat

If all goes well and you plot the result (``gnuplot plot_bonds.gp``), you
should get something like this:

.. image:: Bonds.*

Finally, a third related example.  This one shows how the number of bonds at
an intermediate plate separation varies with solution hybridisation free
energy::

  import dnacc
  from dnacc.units import nm
  import numpy as np
  
  num_tethers = 200
  
  box_L = 200 * nm
  plates = dnacc.Plates(Lx=box_L, Ly=box_L, periodic=True)
  
  for (x, y) in np.random.random_sample((num_tethers/2, 2)):
        plates.add_tether(plate='lower',
                          sticky_end='alpha',
                          L=20 * nm,
                          pos=(x * box_L, y * box_L))
  
  for (x, y) in np.random.random_sample((num_tethers/2, 2)):
        plates.add_tether(plate='upper',
                          sticky_end='alphap',
                          L=20 * nm,
                          pos=(x * box_L, y * box_L))
  
  plates.beta_DeltaG0['alpha', 'alphap'] = 0.0
  plates.at(20 * nm)
  
  print("# Delta G_0 (kT)     Bonds Formed / Max possible")
  for dg in np.linspace(-20.0, 0.0, 21):
        plates.beta_DeltaG0['alpha', 'alphap'] = dg
        plates.update(DeltaG0_only=True)
        print dg, (plates.avg_num_bonds / (num_tethers/2))

This example can be found in ``examples/basic/basic_bonds2.py``.  Run it by
executing the following commands::

  cd examples/basic
  python basic_bonds2.py > bonds2.dat

If all goes well and you plot the result (``gnuplot plot_bonds2.gp``), you
should get something like this:

.. image:: Bonds2.*

Quick overview of what remains
------------------------------

The main class in this module is :class:`.Plates`, used to construct a
system of two parallel DNA-coated plates.  For a mean-field approximation,
use :class:`.PlatesMeanField`.  A more generic interface is exposed in the
:class:`.System` class.

The configurational statistics of the bound DNA tethers are encapsulated in
classes in the :mod:`.tether_statistics` module.  Currently, only
:class:`.RodsGraftedOnPlates`, :class:`.RodsGraftedOnPlatesMeanField` and
:class:`.RodsGraftedOnSpheres` are defined.  See the ``ssDNA_tethers`` and
``comparison_to_lce`` examples for guidance on how to set up different
configurational statistics.

Finally, the :func:`.calc_spheres_potential` function can turn a plate-plate
potential into a sphere-sphere potential using the Derjaguin approximation.

For the moment, the best way to learn how to use this Python module is to
look at the examples in the ``examples`` folder, and to browse through the
:doc:`reference` section of this manual.


.. _dnacc-via-simple-dnacc:

Using DNACC via ``simple_dnacc``
++++++++++++++++++++++++++++++++

``simple_dnacc`` can be run directly as an executable using::

    ./simple_dnacc

and is meant to provide an easy way to run specific calculations with our
python library without the need to write a python program "from scratch".

``simple_dnacc`` reads the ``CONTROL`` file in the current directory (``pwd``)
where instructions about which calculations to perform are given. Read one of
the examples in the ``examples/simple_dnacc`` directory for a some sample input
files. The entries in the ``CONTROL`` file should be almost self-explanatory. For most
wrong choices, the program should stop anyway.

Below is a list of parameters that can be given as input in the
``CONTROL`` file

* ``dg[name_of_strand, name_of_strand]`` (SYMMETRIC DICTIONARY, real number)
  
  Solution hybridisation free energy between strands ``name1``, ``name2``, given
  in kT. Only the couple need be specified (either ``name1, name2`` or
  ``name2, name1``, by default they are set equal). Depending on the type of
  calculation performed, this can be the minimum value for ``dg`` for the
  specific pair or it can be the set value for the whole calculation.
  Example::
  
    dg["alpha", "beta"] = -10

* ``box[0]`` (``box[1]``) (``box`` is an ARRAY) size of x and y dimension of
  the cell to apply Periodic Boundary Conditions.  Example::

    box[0] = 100 * nm
    box[1] = 50 * nm

* ``L[name_of_strand]`` (DICTIONARY, real number) length of strand ``name_of_strand``.
  Note that lengths should be given in the proper unit!  Example::

    L["alpha"] = 20 * nm

* ``R[name_of_sphere]`` (DICTIONARY, real number) radius of the sphere ``name_of_sphere``. 
  Note that lengths should be given in the proper unit!
  Example::

    R["sphere1"] = 100 * nm

* ``sphere_centre[name_of_sphere]`` (DICTIONARY, tuple of 3 real numbers)
  position of the sphere ``name_of_sphere``.  Example::

    sphere_centre["sphere1"] = (10.0 * nm, 0.0 * nm, 0.0 * nm)

* ``np.random.seed(x)`` (numpy function, real) ``x`` is the seed for random numbers

* ``sigma`` (DICTIONARY). ``sigma[plate1][a]`` Density of ``a`` strands on construct ``plate1``

* ``num_type[geometrical_object][type_of_strand]`` (DICTIONARY),
  gives the number of strands on ``geometrical_object`` of type ``type_of_strand``.
  ``geometrical_object`` is an object (sphere or plate) for which the calculation
  is being carried out.  Example::

    num_type["plate1"]["alpha"] = 120

* ``max_num_samples`` (INTEGER VARIABLE), gives the maximum number of
  sampling points when performing a temperature scan. For example, if
  ``calculation`` is set to ``"number of bonds"``, it calculates the number
  of bonds between ``dg = +10`` to ``min(dg)`` at ``max_num_samples`` intervals.
  Example::
  
    max_num_samples=40

* ``explicit`` (LOGICAL VARIABLE), define whether a calculation with explicit
  tethers is performed. Mean-field calculations are performed otherwise.
  Example::
  
    explicit = True

* ``generate_explicit_tethers`` (LOGICAL VARIABLE), define whether a series of
  tethers' positions have to be generated or if they have to be read from a file.
  Example::

    generate_explicit_tethers = True

* ``explicit_tethers_file`` (STRING VARIABLE) name of the file from which the
  positions of the tethers have to be read. If that file is present in the
  current directory, the variable generate_explicit_tethers is set to False
  regardless of was value it is given inside the ``CONTROL`` file.  Example::

    explicit_tethers_file = "TETHERS_CONFIG_FILE"

* ``calculation`` (STRING VARIABLE) type of calculation to be performed. At
  present, possible calculations are

  - ``"sample potential"``: calculate potential at different temperatures up to
    ``min(dg)``
  - ``"number of bonds"``: calculate the number of bonds as a function of ``dg_min``, from
    ``min(dG) = 10`` up to the set value of ``min(dG)``

* ``geometry`` (STRING VARIABLE). Define the geometrical object for the calculation. 
  At the moment, it can be only be ``"spheres"`` or ``"plates"``.  Example::

    geometry = "spheres"

* ``construct`` (STRING VARIABLE). Define the type of tether to be used in the
  calculation. It can either be ``"rods"`` or ``"polymers"``. In the latter case, for the
  time being this simulate a Freely-Jointed-Chain model and can only be used
  together with ``geometry="plates"``.  Example::

    construct = "rods"

* ``output_file`` (STRING VARIABLE). Name of the output file where
  calculation results are stored.  Example::
  
    output_file = "OUTPUT"

