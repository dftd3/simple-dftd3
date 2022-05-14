Comparison with other DFT-D3 implementations
============================================

This DFT-D3 reimplementation was created as a spin-off from the `dftd4`_ and `xtb`_ project, to provide an easier to use API, improve the parallel performance and get a fast build of the DFT-D3 project.
It is however not the only project providing an implementation of DFT-D3, many forks of the original reference implementation and some reimplementations are currently available.

.. _dftd4: https://github.com/dftd4/dftd4
.. _xtb: https://github.com/grimme-lab/xtb

A non-comprehensive list of DFT-D3 implementations is provided here:

============================== ========== ==================== ==========================
 repository                     license    APIs                 notes
============================== ========== ==================== ==========================
 `dftd3`_                       GPL-1.0    Fortran              reference implementation
 `dftd3/simple-dftd3`_          LGPL-3.0   Fortran, C, Python
 `dftbplus/dftd3-lib`_          GPL-1.0    Fortran              patched fork
 `ehermes/ased3`_               LGPL-3.0   Python               f2py, ASE
 `pfnet-research/torch-dftd`_   MIT        Python               torch
 `cuanto/libdftd3`_             GPL-3.0    Fortran, Python      ctypes, pyscf
 `cresset-group/dftd3`_         GPL-1.0    Fortran              patched fork
 `loriab/dftd3`_                GPL-1.0    Fortran              patched fork, Windows
 `f3rmion/dftd3`_               GPL-1.0    Fortran              patched fork
 `bobbypaton/pydftd3`_          MIT        Python               Gaussian
============================== ========== ==================== ==========================

.. _dftd3: http://mctc.uni-bonn.de/software/dft-d3
.. _dftd3/simple-dftd3: https://github.com/dftd3/simple-dftd3
.. _dftbplus/dftd3-lib: https://github.com/dftbplus/dftd3-lib
.. _ehermes/ased3: https://github.com/ehermes/ased3
.. _pfnet-research/torch-dftd: https://github.com/pfnet-research/torch-dftd
.. _cuanto/libdftd3: https://github.com/cuanto/libdftd3
.. _cresset-group/dftd3: https://github.com/cresset-group/dftd3
.. _loriab/dftd3: https://github.com/loriab/dftd3
.. _f3rmion/dftd3: https://github.com/f3rmion/dftd3
.. _bobbypaton/pydftd3: https://github.com/bobbypaton/pyDFTD3

Many more versions are probably around or redistributed in various quantum chemistry programs.


Users of this library
---------------------

A list of projects currently using this DFT-D3 implementation is given here.

`tblite <https://github.com/tblite/tblite>`_: (since 0.1.0)
  Light-weight tight-binding framework
`DFTB+ <https://github.com/dftbplus/dftbplus>`_: (since 21.2)
  General package for performing fast atomistic calculations
`DFT-FE <https://github.com/dftfeDevelopers/dftfe>`_:
  Real-space DFT calculations using Finite Elements
`QCEngine <https://github.com/molssi/qcengine>`_: (WIP)
  Quantum chemistry program executor and IO standardizer.
  For current status see `qcegine#343 <https://github.com/MolSSI/QCEngine/pull/343>`_
`Siesta <https://gitlab.com/siesta-project/siesta>`_: (WIP)
  A first-principles materials simulation code using DFT.
  For current status see `siesta!70 <https://gitlab.com/siesta-project/siesta/-/merge_requests/70>`_


If your project is using *s-dftd3* feel free to add your project to this list.
