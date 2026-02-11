Using 3c composite methods
==========================

In this tutorial we will learn how to use the 3c family of composite electronic structure methods with the ``dftd3`` Python API.
Composite methods like B97-3c\ :footcite:`brandenburg2018` and PBEh-3c\ :footcite:`grimme2015` combine a density functional with a specific basis set, a D3 dispersion correction, and a geometric counter-poise (gCP) correction\ :footcite:`kruse2012` to correct for basis set superposition error (BSSE).


Background on 3c methods
------------------------

The 3c methods are designed as efficient and accurate composite electronic structure methods.
The name "3c" refers to the three corrections applied on top of a base density functional:

1. **D3 dispersion correction** – Accounts for London dispersion interactions using the DFT-D3 model with specifically fitted damping parameters.
2. **Geometric counter-poise (gCP) correction** – Corrects the basis set superposition error that arises from using small basis sets.
3. **Short-range basis set (SRB) correction** – In some methods an additional short-range bond length correction is applied to improve covalent bond lengths.

The total energy of a 3c method is computed as

.. math::

   E_\text{3c} = E_\text{DFT/basis} + E_\text{D3} + E_\text{gCP}

where :math:`E_\text{DFT/basis}` is the SCF energy computed with the corresponding density functional and basis set, :math:`E_\text{D3}` is the DFT-D3 dispersion energy, and :math:`E_\text{gCP}` is the geometric counter-poise correction energy.


Computing the D3 contribution
------------------------------

The D3 dispersion correction is the first component we need.
Each 3c method uses a specific damping scheme and parameters.
For example, PBEh-3c uses rational (BJ) damping with ``s8=0.0``, while B97-3c uses rational damping with its own set of parameters.

Using the Python API we can easily obtain the D3 energy for a given molecule.
In this example we will use a caffeine molecule:

.. code-block:: python

   import numpy as np
   from dftd3.interface import DispersionModel, RationalDampingParam

   # Caffeine molecule (atomic numbers and positions in Bohr)
   numbers = np.array(
       [6, 7, 6, 7, 6, 6, 6, 8, 7, 6, 8, 7, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   )
   positions = np.array([
       [+2.02799738646442, +0.09231312124713, -0.14310895950963],
       [+4.75011007621000, +0.02373496014051, -0.14324124033844],
       [+6.33434307654413, +2.07098865582721, -0.14235306905930],
       [+8.72860718071825, +1.38002919517619, -0.14265542523943],
       [+8.65318821103610, -1.19324866489847, -0.14231527453678],
       [+6.23857175648671, -2.08353643730276, -0.14218299370797],
       [+5.63266886875962, -4.69950321056008, -0.13940509630299],
       [+3.44931709749015, -5.48092386085491, -0.14318454855466],
       [+7.77508917214943, -5.24064112783473, -0.13206210149840],
       [+9.68504443463980, -3.43480556543577, -0.13376503914194],
       [+11.7766597740572, -2.85589272667498, -0.13347836327959],
       [+0.71292821912498, -1.81184541295565, -0.14404507712755],
       [-1.07804988915028, -0.36933811262178, -0.14399838668498],
       [+9.84554065797340, -6.86700842661498, -0.13277505395063],
       [-2.48328028863736, -1.73067674389689, -0.14259442617502],
       [-1.34385752710948, +0.57786478045678, +1.58564153988000],
       [-1.34425478091498, +0.57834894571524, -1.87640258713638],
       [+8.96081504094682, -8.37942090821983, +1.00589803206426],
       [+11.5017680845878, -6.28787412600376, -0.13181456387625],
       [+8.96116647088646, -8.37865627213056, -1.27402975539366],
       [+1.52389863680920, +4.42491864930852, +1.59083729873498],
       [+1.52355529826621, +4.42437909527069, -1.87861372682396],
       [-0.43240687424107, +5.48666820368743, -0.14223727178700],
       [+3.59762505724975, +5.19469189498972, -0.14192321862387],
   ])

   # PBEh-3c D3 parameters (rational damping)
   model = DispersionModel(numbers, positions)
   res = model.get_dispersion(RationalDampingParam(method="pbeh3c"), grad=False)
   print(f"D3 dispersion energy: {res['energy']:16.10f} Hartree")

Since the D3 parameters are already stored in the library for the 3c methods, we can simply pass the method name (e.g. ``"pbeh3c"`` or ``"b973c"``) to ``RationalDampingParam``.

.. note::

   B97-3c also has parameters for zero damping.
   When using it as part of the 3c composite method, the rational damping (BJ) parameters should be used.


Computing the gCP contribution
-------------------------------

The geometric counter-poise (gCP) correction accounts for the basis set superposition error (BSSE) that occurs with small or medium-sized basis sets.
Using the Python API, we can compute the gCP energy using the ``GeometricCounterpoise`` class.

.. code-block:: python

   from dftd3.interface import GeometricCounterpoise

   # Create gCP object for a 3c method
   # For B97-3c, the method name alone is sufficient
   gcp = GeometricCounterpoise(
       numbers,
       positions,
       method="b973c",
   )

   res = gcp.get_counterpoise(grad=False)
   print(f"gCP energy: {res['energy']:16.10f} Hartree")

For methods like PBEh-3c which use a named basis set, both the method and basis set name can be specified:

.. code-block:: python

   gcp = GeometricCounterpoise(
       numbers,
       positions,
       method="pbeh3c",
   )

   res = gcp.get_counterpoise(grad=False)
   print(f"gCP energy: {res['energy']:16.10f} Hartree")

The ``GeometricCounterpoise`` class loads the appropriate gCP parameters internally, including the short-range bond (SRB) correction where applicable.


Putting it all together
-----------------------

To assemble the full 3c composite energy, we combine both the D3 dispersion correction and the gCP correction with the SCF energy from a DFT calculation.

Here is a complete example computing both correction terms for the B97-3c method:

.. code-block:: python

   import numpy as np
   from dftd3.interface import (
       DispersionModel,
       RationalDampingParam,
       GeometricCounterpoise,
   )

   numbers = np.array(
       [6, 7, 6, 7, 6, 6, 6, 8, 7, 6, 8, 7, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   )
   positions = np.array([
       [+2.02799738646442, +0.09231312124713, -0.14310895950963],
       [+4.75011007621000, +0.02373496014051, -0.14324124033844],
       [+6.33434307654413, +2.07098865582721, -0.14235306905930],
       [+8.72860718071825, +1.38002919517619, -0.14265542523943],
       [+8.65318821103610, -1.19324866489847, -0.14231527453678],
       [+6.23857175648671, -2.08353643730276, -0.14218299370797],
       [+5.63266886875962, -4.69950321056008, -0.13940509630299],
       [+3.44931709749015, -5.48092386085491, -0.14318454855466],
       [+7.77508917214943, -5.24064112783473, -0.13206210149840],
       [+9.68504443463980, -3.43480556543577, -0.13376503914194],
       [+11.7766597740572, -2.85589272667498, -0.13347836327959],
       [+0.71292821912498, -1.81184541295565, -0.14404507712755],
       [-1.07804988915028, -0.36933811262178, -0.14399838668498],
       [+9.84554065797340, -6.86700842661498, -0.13277505395063],
       [-2.48328028863736, -1.73067674389689, -0.14259442617502],
       [-1.34385752710948, +0.57786478045678, +1.58564153988000],
       [-1.34425478091498, +0.57834894571524, -1.87640258713638],
       [+8.96081504094682, -8.37942090821983, +1.00589803206426],
       [+11.5017680845878, -6.28787412600376, -0.13181456387625],
       [+8.96116647088646, -8.37865627213056, -1.27402975539366],
       [+1.52389863680920, +4.42491864930852, +1.59083729873498],
       [+1.52355529826621, +4.42437909527069, -1.87861372682396],
       [-0.43240687424107, +5.48666820368743, -0.14223727178700],
       [+3.59762505724975, +5.19469189498972, -0.14192321862387],
   ])

   # 1. Compute D3 dispersion correction
   model = DispersionModel(numbers, positions)
   d3_res = model.get_dispersion(RationalDampingParam(method="b973c"), grad=False)
   d3_energy = d3_res["energy"]

   # 2. Compute gCP correction
   gcp = GeometricCounterpoise(numbers, positions, method="b973c")
   gcp_res = gcp.get_counterpoise(grad=False)
   gcp_energy = gcp_res["energy"]

   print(f"D3 dispersion energy:  {d3_energy:16.10f} Hartree")
   print(f"gCP correction energy: {gcp_energy:16.10f} Hartree")
   print(f"Total correction:      {d3_energy + gcp_energy:16.10f} Hartree")

.. note::

   The SCF energy :math:`E_\text{DFT/basis}` must be computed separately using a quantum chemistry program (such as PySCF, ORCA, or Turbomole) with the correct functional and basis set.
   The total 3c composite energy is then:

   :math:`E_\text{3c} = E_\text{SCF} + E_\text{D3} + E_\text{gCP}`


Including gradients
-------------------

For geometry optimizations, both the D3 and gCP gradients need to be included.
To request gradients, pass ``grad=True``:

.. code-block:: python

   # D3 gradient
   d3_res = model.get_dispersion(RationalDampingParam(method="b973c"), grad=True)
   d3_gradient = d3_res["gradient"]

   # gCP gradient
   gcp_res = gcp.get_counterpoise(grad=True)
   gcp_gradient = gcp_res["gradient"]

   # Total correction gradient
   total_gradient = d3_gradient + gcp_gradient

The gradient arrays have the shape ``(natoms, 3)`` and are in atomic units (Hartree/Bohr).


Summary
-------

In this tutorial we learned how to

- compute the D3 dispersion energy for 3c composite methods using stored parameters
- compute the geometric counter-poise (gCP) correction for basis set superposition error
- combine both corrections to obtain the total correction energy and gradient for 3c methods

The available 3c methods with built-in parameter support include B97-3c, PBEh-3c, and HF-3c.


Literature
----------

.. footbibliography::
