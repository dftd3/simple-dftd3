How to distribute a calculation over MPI ranks?
===============================================

When every MPI rank holds the same structure, each rank would normally repeat the complete dispersion calculation.
A work partition assigns a disjoint share of the pairwise, three-body and reciprocal space loops to each rank instead.

Parts are zero based and every unit of work belongs to exactly one part, so summing the energy, gradient and virial over all parts reproduces the complete result.
D3 itself performs no communication, the reduction is left to the caller.

.. note::

   The :math:`C_6` coefficients are only evaluated for the pairs a part consumes, unless the three-body term is active, which reads them for all pairs.
   The coordination number is consumed in full by every part and stays unpartitioned unless the parts can exchange it, see :ref:`reducer`.

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


Let D3 do the reduction
-----------------------

MPI support is an opt-in build option:

.. tab-set::

   .. tab-item:: meson

      .. code-block:: shell

         meson setup _build -Dmpi=true

   .. tab-item:: CMake

      .. code-block:: shell

         cmake -B _build -DSDFTD3_WITH_MPI=ON

It enables the ``dftd3_mpi`` module, which derives the partition from a communicator and reduces the results, replacing the whole boilerplate above.
The C API takes the communicator on the model instead of a part index, everything else stays the same:

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use dftd3, only : get_dispersion_mpi
         use mpi, only : MPI_COMM_WORLD

         call get_dispersion_mpi(error, mol, disp, param, realspace_cutoff(), &
            & MPI_COMM_WORLD, energy, gradient)

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         dftd3_set_model_mpi_comm(error, disp, MPI_Comm_c2f(MPI_COMM_WORLD));
         dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, NULL);

Communicators are passed as plain Fortran handles, users of ``mpi_f08`` pass ``comm%mpi_val`` and C callers convert with ``MPI_Comm_c2f``.
The counter-poise correction has the same entry points, ``get_counterpoise_mpi`` and ``dftd3_set_gcp_mpi_comm``.
For everything else, ``new_mpi_work_partition`` creates the partition of the calling rank and leaves the reduction to you.

.. note::

   The routines are collective, every rank of the communicator has to call them with the same input.

Without MPI support the entry points are still there, but every one of them reports the missing feature in the error handler instead.
The build configuration can be queried at compile time and at runtime:

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use dftd3, only : dftd3_has_mpi, dftd3_has_feature

         if (dftd3_has_feature("mpi")) then
            ! ...
         end if

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         if (dftd3_has_feature("mpi")) {
            /* ... */
         }

The command line driver reports the same in ``s-dftd3 --version``.


.. _reducer:

Partition the coordination number
---------------------------------

The coordination number is needed in full by every part, so it can only be partitioned if the parts exchange it halfway through the calculation.
Passing a ``work_reducer`` alongside the partition enables this, ``get_dispersion_mpi`` does so automatically.
Bring your own communication layer by implementing the deferred ``reduce`` binding:

.. code-block:: fortran

   type, extends(work_reducer) :: my_reducer
   contains
      procedure :: reduce => my_reduce
   end type my_reducer

   ! sum val over all parts, in place
   subroutine my_reduce(self, val, error)
      class(my_reducer), intent(in) :: self
      real(wp), intent(inout), contiguous :: val(:)
      type(error_type), allocatable, intent(out) :: error
   end subroutine my_reduce

The reducer is called twice per gradient evaluation, both times on an array of the size of the system.
Without it the coordination number stages stay unpartitioned and the results are unchanged.
