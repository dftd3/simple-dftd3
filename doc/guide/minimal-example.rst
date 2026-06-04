How to use this library?
========================

This section contains a few self-contained examples on how to use D3.


Compute energy with rational damping
------------------------------------

This example shows how to compute the dispersion energy with the rational damping function.

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. literalinclude:: minimal-example/energy.f90
         :language: fortran
         :caption: energy.f90

   .. tab-item:: C
      :sync: c

      .. literalinclude:: minimal-example/energy.c
         :language: c
         :caption: energy.c

   .. tab-item:: Python
      :sync: python

      .. literalinclude:: minimal-example/energy.py
         :language: python
         :caption: energy.py

To test this example you can install the dependencies with

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

You can run the example code with

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: shell

         ❯ $FC energy.f90 $(pkg-config s-dftd3 mctc-lib --cflags --libs) && ./a.out
         Dispersion energy for PBE0-D3(BJ) is -0.0009218696 Hartree

   .. tab-item:: C
      :sync: c

      .. code-block:: shell

         ❯ $CC energy.c $(pkg-config s-dftd3 mctc-lib --cflags --libs) && ./a.out
         Dispersion energy for PBE0-D3(BJ) is -0.0009218696 Hartree

   .. tab-item:: Python
      :sync: python

      .. code-block:: shell

         ❯ python energy.py
         Dispersion energy for PBE0-D3(BJ) is -0.0009219059 Hartree


Disable selected atoms with ghost indices
-----------------------------------------

The same idea works for all public interfaces. Use the atom indices of the
fragment that should be excluded from the dispersion contribution, then pass
them to the model or CLI entry point:

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. literalinclude:: minimal-example/ghost.f90
         :language: fortran
         :caption: ghost.f90

   .. tab-item:: C
      :sync: c

      .. literalinclude:: minimal-example/ghost.c
         :language: c
         :caption: ghost.c

   .. tab-item:: Python
      :sync: python

      .. literalinclude:: minimal-example/ghost.py
         :language: python
         :caption: ghost.py

The Fortran constructor uses 1-based atom indices, while the C and Python
interfaces use 0-based indices.
For the CLI, the same selection is available as

.. code-block:: shell

   ❯ s-dftd3 run --ghost 4,5,6 structure.xyz
