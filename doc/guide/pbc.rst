How to run periodic calculations?
=================================

For a periodic system the two-body dispersion energy is a lattice sum, which has to be either truncated in real space or transformed into reciprocal space.
This page shows the three available setups on the same system, the cubic phase of ammonia with r²SCAN-D3(BJ).

.. note::

   The background on the truncation error and how to size a cutoff for a target accuracy is covered in :doc:`accuracy`.

To test the examples you can install the dependencies with

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: text

         mamba create d3 simple-dftd3 fortran-compiler pkg-config
         mamba activate d3

   .. tab-item:: C
      :sync: c

      .. code-block:: text

         mamba create d3 simple-dftd3 c-compiler pkg-config
         mamba activate d3

   .. tab-item:: Python
      :sync: python

      .. code-block:: text

         mamba create d3 dftd3-python
         mamba activate d3


Remove the truncation error with reciprocal space summation
-----------------------------------------------------------

Enabling the low-rank expansion of the reference :math:`C_6` coefficients makes the two-body energy separable, so it can be summed over the reciprocal lattice instead.
The damped pair potential is bounded at the origin, which lets the reciprocal sum converge exponentially without a real space complement, and the two-body cutoff drops out of the calculation entirely.

This requires three-dimensional periodic boundary conditions and a damping function with a known reciprocal space representation, currently the rational and the zero damping function.

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. literalinclude:: pbc-example/ewald.f90
         :language: fortran
         :caption: ewald.f90

   .. tab-item:: C
      :sync: c

      .. literalinclude:: pbc-example/ewald.c
         :language: c
         :caption: ewald.c

   .. tab-item:: Python
      :sync: python

      .. literalinclude:: pbc-example/ewald.py
         :language: python
         :caption: ewald.py

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: shell

         ❯ $FC ewald.f90 $(pkg-config s-dftd3 mctc-lib --cflags --libs) && ./a.out
         Dispersion energy for r2SCAN-D3(BJ) is -0.0122603293 Hartree

   .. tab-item:: C
      :sync: c

      .. code-block:: shell

         ❯ $CC ewald.c $(pkg-config s-dftd3 --cflags --libs) && ./a.out
         Dispersion energy for r2SCAN-D3(BJ) is -0.0122603293 Hartree

   .. tab-item:: Python
      :sync: python

      .. code-block:: shell

         ❯ python ewald.py
         Dispersion energy for r2SCAN-D3(BJ) is -0.0122603293 Hartree

The same is available from the command line with

.. code-block:: shell

   ❯ s-dftd3 run --ewald --bj r2scan structure.poscar

Note that the pairwise decomposition and the analytical Hessian are only implemented for the real space summation and are rejected for such a model.


Choose how the reciprocal sum is evaluated
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default the structure factors are interpolated on a mesh with cardinal B-splines and the sum is evaluated with a fast Fourier transform, which costs :math:`O(N \log N)`.
Summing over the reciprocal lattice directly instead visits every atom at every wave vector and costs :math:`O(N^2)`.

The mesh is derived from the reciprocal cutoff unless it is given explicitly:

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         ! explicit mesh, rounded up to a power of two
         call new_d3_model(disp, mol, lowrank=d3_lowrank_config(mesh=32))

         ! sum over the reciprocal lattice directly
         call new_d3_model(disp, mol, lowrank=d3_lowrank_config(mesh=-1))

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         /* explicit mesh, rounded up to a power of two */
         dftd3_set_model_ewald(error, d3, 0, 0.0, 0.0, 32);

         /* sum over the reciprocal lattice directly */
         dftd3_set_model_ewald(error, d3, 0, 0.0, 0.0, -1);

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         # explicit mesh, rounded up to a power of two
         model.set_ewald_summation(mesh=32)

         # sum over the reciprocal lattice directly
         model.set_ewald_summation(mesh=-1)

Both reach the same energy for the cell above

.. code-block:: text

   Dispersion energy for r2SCAN-D3(BJ) is -0.0122603293 Hartree

From the command line the mesh is selected with ``--ewald-mesh`` and the direct summation with ``--ewald-direct``

.. code-block:: shell

   ❯ s-dftd3 run --ewald-mesh 32 --bj r2scan structure.poscar
   ❯ s-dftd3 run --ewald-direct --bj r2scan structure.poscar

See :doc:`accuracy` for how to choose the mesh.


Smooth the real space cutoff
----------------------------

A hard cutoff makes the energy jump whenever an atom pair crosses it, which shows up as noise in a geometry optimization or a lattice scan.
Giving the cutoff a smoothing width switches the pair contributions off over the last fraction of a Bohr instead, so the energy stays continuous.
The width is the trailing argument for the two- and three-body cutoff and only changes the energy by the size of the contributions it damps.

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. literalinclude:: pbc-example/smooth.f90
         :language: fortran
         :caption: smooth.f90

   .. tab-item:: C
      :sync: c

      .. literalinclude:: pbc-example/smooth.c
         :language: c
         :caption: smooth.c

   .. tab-item:: Python
      :sync: python

      .. literalinclude:: pbc-example/smooth.py
         :language: python
         :caption: smooth.py

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: shell

         ❯ $FC smooth.f90 $(pkg-config s-dftd3 mctc-lib --cflags --libs) && ./a.out
         Dispersion energy for r2SCAN-D3(BJ) is -0.0122446362 Hartree

   .. tab-item:: C
      :sync: c

      .. code-block:: shell

         ❯ $CC smooth.c $(pkg-config s-dftd3 --cflags --libs) && ./a.out
         Dispersion energy for r2SCAN-D3(BJ) is -0.0122446362 Hartree

   .. tab-item:: Python
      :sync: python

      .. code-block:: shell

         ❯ python smooth.py
         Dispersion energy for r2SCAN-D3(BJ) is -0.0122446362 Hartree

Smoothing removes the discontinuity, it does not remove the truncation error.
The result still differs from the reciprocal space value by about 1.6 · 10⁻⁵ Hartree for this cell.


Truncate the real space sum
---------------------------

Without any further setup the lattice sum is truncated at the default cutoffs, 60 Bohr for the two-body term and 40 Bohr for the three-body term and the coordination number.
This is the cheapest option and the one used by every non-periodic interface.

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. literalinclude:: pbc-example/energy.f90
         :language: fortran
         :caption: energy.f90

   .. tab-item:: C
      :sync: c

      .. literalinclude:: pbc-example/energy.c
         :language: c
         :caption: energy.c

   .. tab-item:: Python
      :sync: python

      .. literalinclude:: pbc-example/energy.py
         :language: python
         :caption: energy.py

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: shell

         ❯ $FC energy.f90 $(pkg-config s-dftd3 mctc-lib --cflags --libs) && ./a.out
         Dispersion energy for r2SCAN-D3(BJ) is -0.0122446595 Hartree

   .. tab-item:: C
      :sync: c

      .. code-block:: shell

         ❯ $CC energy.c $(pkg-config s-dftd3 --cflags --libs) && ./a.out
         Dispersion energy for r2SCAN-D3(BJ) is -0.0122446595 Hartree

   .. tab-item:: Python
      :sync: python

      .. code-block:: shell

         ❯ python energy.py
         Dispersion energy for r2SCAN-D3(BJ) is -0.0122446595 Hartree

Raising the two-body cutoff reduces the remaining error, see :doc:`accuracy` for choosing a cutoff for a target accuracy.
