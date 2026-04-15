Real-space cutoffs for periodic calculations
============================================

For periodic structures D3 sums pair and, optionally, ATM three-body
contributions over lattice translations in real space. When comparing the
standalone executable, the Python API, the ASE calculator, or another code
such as VASP, check these settings first:

* whether the ATM term is enabled,
* whether the same real-space cutoffs are used,
* and whether the reported energies are compared in the same unit.

The library defaults are ``disp2 = 60 Bohr``, ``disp3 = 40 Bohr``, and
``cn = 40 Bohr``.


Units and defaults
------------------

+----------------------+----------------------------------------+-----------------------+
| Interface            | Real-space cutoff input                | Returned energy       |
+======================+========================================+=======================+
| ``s-dftd3``          | not configurable from the CLI; uses    | Hartree (``Eh``)      |
|                      | the library defaults above             |                       |
+----------------------+----------------------------------------+-----------------------+
| ``dftd3.interface``  | Bohr                                   | Hartree               |
+----------------------+----------------------------------------+-----------------------+
| ``dftd3.ase.DFTD3``  | Angstrom; use ``60 * Bohr`` to pass    | eV, following ASE     |
|                      | 60 Bohr through ASE                    | conventions           |
+----------------------+----------------------------------------+-----------------------+

The standalone executable treats periodic input formats such as ``POSCAR`` or
``CONTCAR`` automatically as periodic structures. VASP reports energies in eV,
so convert the ``s-dftd3`` or Python API energies before comparing them.


Standalone executable
---------------------

The executable uses the default periodic cutoffs. Without ``--atm`` it computes
the two-body D3 energy; adding ``--atm`` enables the ATM term.

.. code-block:: shell

   s-dftd3 CONTCAR --bj pbe
   s-dftd3 CONTCAR --bj pbe --atm

If you need to study the cutoff dependence explicitly, use the Python API or
the ASE calculator below. Those interfaces expose the periodic cutoffs
directly.


Python API
----------

The low-level Python API expects coordinates, lattice vectors, and cutoffs in
Bohr and returns the energy in Hartree.

.. code-block:: python

   from ase.io import read
   from ase.units import Bohr, Hartree
   from dftd3.interface import DispersionModel, RationalDampingParam

   atoms = read("CONTCAR")
   model = DispersionModel(
       atoms.numbers,
       atoms.positions / Bohr,
       lattice=atoms.cell.array / Bohr,
       periodic=atoms.pbc,
   )

   params = {
       "D3(BJ)": RationalDampingParam(method="PBE", atm=False),
       "D3(BJ)-ATM": RationalDampingParam(method="PBE", atm=True),
   }

   for disp2, disp3, cn in [(60.0, 40.0, 40.0), (95.0, 70.0, 70.0)]:
       model.set_realspace_cutoff(disp2=disp2, disp3=disp3, cn=cn)
       for label, param in params.items():
           energy = model.get_dispersion(param=param, grad=False)["energy"]
           print(
               label,
               f"cutoff=({disp2}, {disp3}, {cn}) Bohr",
               f"{energy:.12f} Eh",
               f"{energy * Hartree:.8f} eV",
           )

Using ``atm=False`` reproduces the executable call without ``--atm``. The
default ``RationalDampingParam(method=\"PBE\")`` uses the ATM-enabled
parameterization.


ASE calculator
--------------

The ASE calculator expects cutoff values in Angstrom and returns energies in eV.
To pass the same numerical cutoff as the library defaults, multiply by
``ase.units.Bohr``.

.. code-block:: python

   from ase.io import read
   from ase.units import Bohr
   from dftd3.ase import DFTD3

   atoms = read("CONTCAR")

   setups = {
       "D3(BJ)": {"params_tweaks": {"method": "PBE", "atm": False}},
       "D3(BJ)-ATM": {"params_tweaks": {"method": "PBE", "atm": True}},
   }

   cutoffs = [
       {"disp2": 60 * Bohr, "disp3": 40 * Bohr, "cn": 40 * Bohr},
       {"disp2": 95 * Bohr, "disp3": 70 * Bohr, "cn": 70 * Bohr},
   ]

   for label, kwargs in setups.items():
       for cutoff in cutoffs:
           atoms.calc = DFTD3(
               damping="d3bj",
               realspace_cutoff=cutoff,
               cache_api=False,
               **kwargs,
           )
           print(label, cutoff, f"{atoms.get_potential_energy():.8f} eV")

``DFTD3(method=\"PBE\", damping=\"d3bj\")`` uses ATM-enabled parameters. To
match ``s-dftd3 CONTCAR --bj pbe``, request the pairwise-only parameterization
with ``params_tweaks={\"method\": \"PBE\", \"atm\": False}``.
If you prefer explicit damping parameters, setting ``s9 = 0.0`` in
``params_tweaks`` disables the ATM contribution as well.
