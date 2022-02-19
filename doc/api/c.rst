C API
=====

The C API bindings are provided by using the ``iso_c_binding`` intrinsic module.
Generally, objects are exported as opaque pointers and can only be manipulated within the library.
The API user is required delete all objects created in the library by using the provided deconstructor functions to avoid mamory leaks.

Overall four classes of objects are provided by the library

- error handlers (``dftd3_error``),
  used to communicate exceptional conditions and errors from the library to the user
- structure containers (``dftd3_structure``),
  used to represent the system specific information and geometry data,
  only the latter are mutable for the user
- dispersion model objects (``dftd3_model``),
  general model for calculating dispersion releated properties
- damping function objects (``dftd3_param``)
  polymorphic objects to represent the actual method parametrisation

.. note::

   Generally, all quantities provided to the library are assumed to be in `atomic units <https://en.wikipedia.org/wiki/Hartree_atomic_units>`_.

.. contents::


Error handling
--------------

The library provides a light error handle type (``dftd3_error``) for storing error information
The error handle requires only small overhead to construct and can only contain a single error.

The handler is represented by an opaque pointer and can only be manipulated by call from the library.
The user of those objects is required to delete the handlers again using the library provided deconstructors to avoid memory leaks.


Structure data
--------------

The structure data is used to represent the system of interest in the library.
It contains immutable system specific information like the number of atoms, the unique atom groups and the boundary conditions as well as mutable geometry data like cartesian coordinates and lattice parameters.


Performing calculations
-----------------------

An example wrapper to perform a DFT-D3(BJ)-ATM calculation is shown below.


.. code-block:: c

   #include <stdbool.h>
   #include <stdio.h>
   #include <stdlib.h>

   #include "dftd3.h"

   static const buffersize = 512;

   int
   calc_dftd3(int natoms, int* numbers, double* positions,
              double* lattice, bool* periodic, char* method,
              double* energy, double* gradient, double* sigma)
   {
     // Local API objects from the s-dftd3 library
     dftd3_error error = dftd3_new_error();
     dftd3_structure mol = NULL;
     dftd3_model disp = NULL;
     dftd3_param param = NULL;
     int stat = EXIT_SUCCESS;

     // Create a new geometry for the library to work with
     mol = dftd3_new_structure(error, natoms, numbers, positions, lattice, periodic);
     stat = dftd3_check_error(error);

     if (stat) {
       // Initialize the D3 dispersion model for the given structure,
       // this step depends on the atomic numbers, but not on the actual geometry
       disp = dftd3_new_d3_model(error, mol);
       stat = dftd3_check_error(error);
     }

     if (stat) {
       // Load D3(BJ)-ATM parameters for the given method from internal storage,
       // this step depends on the atomic numbers, but not on the actual geometry
       param = dftd3_load_rational_damping(error, mol, method, true);
       stat = dftd3_check_error(error);
     }

     if (stat) {
       // Evaluate the dispersion energy, gradient and virial,
       // the gradient and virial are optional and can be replaced by NULL
       dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
       stat = dftd3_check_error(error);
     }

     if (!stat) {
       char buffer[buffersize];
       dftd3_get_error(error, buffer, buffersize);
       printf("[Error] %s\n", buffer);
     }

     // Always free the used memory
     dftd3_delete(error);
     dftd3_delete(mol);
     dftd3_delete(disp);
     dftd3_delete(param);

     return stat;
   }
