Periodic systems: cutoffs and summation accuracy
================================================

For a periodic system the two-body dispersion energy is a lattice sum.
Because the leading term decays only as :math:`r^{-6}`, truncating the sum at a real space cutoff leaves a systematic error which converges slowly with the cutoff radius.
This page documents the size of that error, how to pick a cutoff for a target accuracy, and how to remove it entirely with the reciprocal space summation.


Truncation error of the real space sum
--------------------------------------

Beyond the damping region the pair interaction approaches the bare :math:`C_6` tail.
Integrating it over the :math:`d` periodic dimensions gives the energy missing from a sum truncated at :math:`R`

.. math::

   \frac{\Delta E}{N} = \frac{\rho\,\langle C_6\rangle\, S_d}{2\,(6-d)\,R^{6-d}}
   \qquad
   \langle C_6\rangle = \frac{1}{N^2}\sum_{ij} C_6^{ij}

with the number density :math:`\rho` of the periodic cell and the surface :math:`S_d` of the :math:`d`-dimensional unit sphere.
For a bulk system this is :math:`\Delta E/N = 2\pi\rho\langle C_6\rangle/(3R^3)`, i.e. the error decays only as :math:`R^{-3}`.
The relation is independent of the damping function, since every damping function approaches one at large distances.

The table below lists the error per atom for PBE-D3(BJ) against a converged reciprocal space reference.
All values in Hartree per atom.

============ ============ ============ ============ ============ ============
 disp2 / a₀   benzene      urea         anthracene   ice Ih       ice VII
============ ============ ============ ============ ============ ============
 20           4.6e-5       3.7e-5       5.6e-5       2.0e-5       3.0e-5
 30           1.4e-5       1.1e-5       1.6e-5       5.9e-6       9.1e-6
 40           5.8e-6       4.6e-6       6.9e-6       2.5e-6       3.8e-6
 60           1.7e-6       1.4e-6       2.0e-6       7.3e-7       1.1e-6
 80           7.2e-7       5.8e-7       8.6e-7       3.1e-7       4.8e-7
 100          3.7e-7       3.0e-7       4.4e-7       1.6e-7       2.4e-7
 140          1.3e-7       1.1e-7       1.6e-7       5.8e-8       8.9e-8
 200          4.6e-8       3.7e-8       5.5e-8       2.0e-8       3.0e-8
 300          1.4e-8       1.1e-8       1.6e-8       5.9e-9       9.0e-9
============ ============ ============ ============ ============ ============

The measured amplitudes :math:`\Delta E R^3/N` agree with the predicted :math:`2\pi\rho\langle C_6\rangle/3` to better than one percent:

=============== ============ ============ ============ ============ ============
                 benzene      urea         anthracene   ice Ih       ice VII
=============== ============ ============ ============ ============ ============
 predicted       0.369        0.295        0.439        0.158        0.244
 measured        0.369        0.295        0.439        0.158        0.244
=============== ============ ============ ============ ============ ============

The default cutoff of 60 a₀ therefore leaves an error of roughly 1 to 2 µHartree per atom for a molecular crystal, which is small on an absolute scale but can be comparable to the energy differences between polymorphs.
Note that the cost of the real space sum grows as :math:`R^3`, so reducing the error by one order of magnitude costs about a factor of ten in time.


Choosing a cutoff for a target accuracy
---------------------------------------

Inverting the relation above yields the cutoff needed for a requested accuracy.
The ``get_realspace_cutoff`` function performs this estimate for a given structure and dispersion model:

.. code-block:: fortran

   use dftd3, only : realspace_cutoff, get_realspace_cutoff

   type(realspace_cutoff) :: cutoff

   ! two-body energy converged to 1 µHartree per atom
   cutoff = get_realspace_cutoff(mol, disp, 1.0e-6_wp)

The estimate covers three, two and one dimensional boundary conditions, and returns a cutoff spanning the whole system for a finite one.
It assumes a scaling factor :math:`s_6` of one and neglects the faster decaying :math:`C_8` tail, both of which are covered by a safety margin.

Only the two-body cutoff is derived from the accuracy target.
The three-body cutoff is left at its default, and the coordination number cutoff is not an accuracy parameter at all: the D3 counting function has a non-vanishing limit at large distances, so increasing ``cn`` does not converge the coordination number but changes it.
Energies computed with different ``cn`` values are not comparable at the sub-µHartree level, which matters when validating cutoff settings against each other.


Reciprocal space summation
--------------------------

The truncation error can be removed altogether by evaluating the two-body sum in reciprocal space.\ :footcite:`valeeva2026`
This requires a separable representation of the environment dependent :math:`C_6` coefficients, which is created by passing a ``d3_lowrank_config`` to the model constructor:

.. code-block:: fortran

   use dftd3, only : d3_model, new_d3_model, d3_lowrank_config

   type(d3_model) :: disp

   call new_d3_model(disp, mol, lowrank=d3_lowrank_config())

A model set up this way uses the Ewald summation for the two-body energy whenever the system is periodic in all three dimensions and the damping function supports it, and ``disp2`` is then ignored.
The rational and zero damping functions are supported; all others report an error.

The accuracy of the expansion and of the reciprocal sum are controlled independently:

``tolerance``
  Maximum relative error of the reconstructed reference :math:`C_6` coefficients.
  The default of 1e-4 needs a rank below ten for typical systems.

``rank``
  Fixed rank of the expansion, overriding ``tolerance``.

``kcut``
  Reciprocal space cutoff in inverse Bohr.
  The default of zero derives it from the damping radii.

``mesh``
  Number of particle mesh points per direction, see below.
  The default of zero derives the mesh from ``kcut``, a negative value sums over the reciprocal lattice directly.

Convergence with respect to ``kcut`` is exponential.
The table lists the deviation from a reference computed with ``kcut`` of 14 a₀\ :sup:`-1`, in Hartree per atom:

============ ============ ============ ============ ============ ============
 kcut         benzene      urea         anthracene   ice Ih       ice VII
============ ============ ============ ============ ============ ============
 4            2.7e-09      1.4e-08      3.1e-10      1.1e-08      1.9e-08
 6            6.5e-11      2.0e-11      7.8e-11      7.9e-12      5.4e-11
 8            7.8e-13      1.6e-13      6.4e-13      3.9e-13      3.6e-13
 10           1.9e-14      2.1e-14      1.7e-14      1.6e-14      2.6e-14
============ ============ ============ ============ ============ ============

A cutoff of 8 a₀\ :sup:`-1` reaches one picohartree per atom, and the automatic setting is more accurate than any real space cutoff of practical size.
For the zero damping function the transform decays more slowly and roughly twice the cutoff is needed for the same accuracy.

The converged two-body energies used as reference above, for PBE-D3(BJ) with :math:`s_9 = 0` and the default coordination number cutoff, are

=================== ======== ==========================
 structure           atoms    E(2) / Hartree
=================== ======== ==========================
 X23 benzene         48       -1.26322976499534e-01
 X23 urea            16       -3.53548720665777e-02
 X23 anthracene      48       -1.49958203449099e-01
 ICE10 ice Ih        48       -5.87846231499417e-02
 ICE10 ice VII       48       -9.93223913624172e-02
=================== ======== ==========================

Reproducing these with the real space summation requires a cutoff beyond 300 a₀.


Particle mesh evaluation
------------------------

The direct reciprocal space sum visits every atom at every wave vector, so its cost grows quadratically with the number of atoms.
Setting ``mesh`` interpolates the structure factors on an equispaced mesh with cardinal B-splines and evaluates the sum with a fast Fourier transform instead, which scales as :math:`O(N \log N)`.
This is the default; a negative ``mesh`` falls back to the direct summation.

The discretisation error is governed by the mesh spacing, not by the number of mesh points, so the mesh has to grow with the cell.
The table lists the deviation from the direct summation for cubic carbon monoxide supercells, in Hartree per atom, against the resulting spacing in a₀:

============ ============ ============
 spacing      atoms        deviation
============ ============ ============
 1.25         250          8.1e-07
 1.00         128          2.0e-07
 0.75         54           5.7e-09
 0.63         250          2.1e-10
 0.50         128          2.4e-10
 0.38         54           2.8e-11
============ ============ ============

A spacing near one Bohr matches a well converged real space calculation, and half a Bohr brings the mesh error below the reciprocal cutoff error.

The two settings interact, because the mesh spans the whole Brillouin zone of the mesh while the direct sum stops at ``kcut``.
Refining the mesh therefore converges past the direct summation rather than towards it.

The cost is set by the mesh rather than by the number of atoms, so it is nearly flat in system size.
Wall times for the two-body energy of the same supercells, on four threads:

============ ============ ============ ============ ============
 atoms        real space   direct       mesh 32³     mesh 64³
============ ============ ============ ============ ============
 54           0.002        0.126        0.045        0.346
 128          0.005        0.565        0.069        0.354
 250          0.010        1.966        0.048        0.381
============ ============ ============ ============ ============

Against the direct summation the mesh is a clear win and grows more so with system size.
Against the truncated real space sum it is not: that sum is :math:`O(N)` with a small prefactor and stays cheaper at the accuracy its default cutoff delivers.
The mesh earns its cost where the real space cutoff cannot reach the required accuracy, which for the systems in the table above means anything below roughly 1e-6 Hartree per atom.

Under a work partition the terms of the low-rank expansion are handed out whole, so the parallelism of the mesh evaluation is capped by the rank rather than by the size of the system.
For a three-element crystal that rank is 5 at the default tolerance and 9 at a tolerance of 1e-6, and parts beyond it stay idle.
The direct summation partitions over wave vectors instead and keeps scaling well past that, which is a reason to prefer it when distributing a small cell over many ranks.

.. footbibliography::
