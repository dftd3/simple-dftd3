How to distribute a calculation over MPI ranks?
===============================================

When every MPI rank holds the same structure, each rank would normally repeat the complete dispersion calculation.
A work partition assigns a disjoint share of the pairwise, three-body and reciprocal space loops to each rank instead.

Parts are zero based and every unit of work belongs to exactly one part, so summing the energy, gradient and virial over all parts reproduces the complete result.
D3 itself performs no communication, the reduction is left to the caller.

.. note::

   Structure dependent quantities such as coordination numbers and :math:`C_6` coefficients are evaluated for the full system on every part.
   The speedup is therefore bound by the interaction loops, which dominate for larger systems.

To test the examples you can install the dependencies with

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: text

         mamba create d3 simple-dftd3 fortran-compiler openmpi pkg-config
         mamba activate d3

   .. tab-item:: C
      :sync: c

      .. code-block:: text

         mamba create d3 simple-dftd3 c-compiler openmpi pkg-config
         mamba activate d3


Partition the interaction loops
-------------------------------

The Fortran API takes the partition as an optional argument, omitting it selects the complete work.
The C API stores it on the dispersion model, next to the cutoffs.

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. literalinclude:: partition-example/mpi.f90
         :language: fortran
         :caption: mpi.f90

   .. tab-item:: C
      :sync: c

      .. literalinclude:: partition-example/mpi.c
         :language: c
         :caption: mpi.c

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: shell

         ❯ mpifort mpi.f90 $(pkg-config s-dftd3 mctc-lib --cflags --libs) && mpirun -n 4 ./a.out
         Dispersion energy for PBE0-D3(BJ) is -0.0009218696 Hartree

   .. tab-item:: C
      :sync: c

      .. code-block:: shell

         ❯ mpicc mpi.c $(pkg-config s-dftd3 --cflags --libs) && mpirun -n 4 ./a.out
         Dispersion energy for PBE0-D3(BJ) is -0.0009218696 Hartree

The result is independent of the number of ranks up to the summation order.
Running with a single rank is identical to omitting the partition, which is also available as ``serial_work_partition`` if you prefer to always pass it explicitly.

The analytical Hessian and the geometric counterpoise correction accept the same partition, the pairwise decomposition does not.
