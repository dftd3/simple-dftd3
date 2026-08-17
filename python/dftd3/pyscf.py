# This file is part of s-dftd3.
# SPDX-Identifier: LGPL-3.0-or-later
#
# s-dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# s-dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.
"""
PySCF Support
=============

Compatibility layer for supporting DFT-D3 and GCP in `pyscf <https://pyscf.org/>`_.
"""

try:
    from pyscf import gto, lib, pbc, mcscf, scf
    from pyscf.grad import rhf as rhf_grad
except ModuleNotFoundError:
    raise ModuleNotFoundError("This submodule requires pyscf installed")

import numpy as np
from typing import Dict, Optional, Tuple, Union

from .interface import (
    DispersionModel,
    RationalDampingParam,
    ZeroDampingParam,
    ModifiedRationalDampingParam,
    ModifiedZeroDampingParam,
    OptimizedPowerDampingParam,
    CSODampingParam,
    ZDampingParam,
    GeometricCounterpoise,
)

GradientsBase = getattr(rhf_grad, "GradientsBase", rhf_grad.Gradients)

_damping_param = {
    "d3bj": RationalDampingParam,
    "d3zero": ZeroDampingParam,
    "d3bjm": ModifiedRationalDampingParam,
    "d3mbj": ModifiedRationalDampingParam,
    "d3zerom": ModifiedZeroDampingParam,
    "d3mzero": ModifiedZeroDampingParam,
    "d3op": OptimizedPowerDampingParam,
    "d3cso": CSODampingParam,
    "d3z": ZDampingParam,
}


class DFTD3Dispersion(lib.StreamObject):
    """
    Implementation of the interface for using DFT-D3 in pyscf.
    The `xc` functional can be provided in the constructor together with the
    `version` of the DFT-D3 damping function to use.
    Possible damping functions are

    ``"d3bj"``: (default)
        For rational damping function
    ``"d3zero"``
        For zero damping function
    ``"d3mbj"``
        Modified damping parameters for the rational damping function
    ``"d3mzero"``
        Modified version of the zero damping function
    ``"d3op"``
        Optimized power damping function
    ``"d3cso"``
        CSO (C6-scaled only) damping function
    ``"d3z"``
        Z damping function

    Custom parameters can be provided with the `param` dictionary.
    The `param` dict contains the damping parameters, at least s8, a1 and a2
    must be provided for rational damping, while s8 and rs6 are required in case
    of zero damping. For CSO and Z damping, a1 must be provided.

    Parameters for (modified) rational damping are:

    ======================== =========== ============================================
    Tweakable parameter      Default     Description
    ======================== =========== ============================================
    s6                       1.0         Scaling of the dipole-dipole dispersion
    s8                       None        Scaling of the dipole-quadrupole dispersion
    s9                       1.0         Scaling of the three-body dispersion energy
    a1                       None        Scaling of the critical radii
    a2                       None        Offset of the critical radii
    alp                      14.0        Exponent of the zero damping (ATM only)
    ======================== =========== ============================================

    Parameters for (modified) zero damping are:

    ======================== =========== ===================================================
    Tweakable parameter      Default     Description
    ======================== =========== ===================================================
    s6                       1.0         Scaling of the dipole-dipole dispersion
    s8                       None        Scaling of the dipole-quadrupole dispersion
    s9                       1.0         Scaling of the three-body dispersion energy
    rs6                      None        Scaling of the dipole-dipole damping
    rs8                      1.0         Scaling of the dipole-quadrupole damping
    alp                      14.0        Exponent of the zero damping
    bet                      None        Offset for damping radius (modified zero damping)
    ======================== =========== ===================================================

    Parameters for optimized power damping are:

    ======================== =========== ============================================
    Tweakable parameter      Default     Description
    ======================== =========== ============================================
    s6                       1.0         Scaling of the dipole-dipole dispersion
    s8                       None        Scaling of the dipole-quadrupole dispersion
    s9                       1.0         Scaling of the three-body dispersion energy
    a1                       None        Scaling of the critical radii
    a2                       None        Offset of the critical radii
    alp                      14.0        Exponent of the zero damping (ATM only)
    bet                      None        Power for the zero-damping component
    ======================== =========== ============================================

    Parameters for CSO (C6-scaled only) damping are:

    ======================== =========== ============================================
    Tweakable parameter      Default     Description
    ======================== =========== ============================================
    s6                       1.0         Scaling of the dipole-dipole dispersion
    s9                       1.0         Scaling of the three-body dispersion energy
    a1                       None        Sigmoid amplitude parameter
    a2                       2.5         Sigmoid reference distance scale
    a3                       0.0         Denominator critical radii scale
    a4                       6.25        Denominator constant offset
    alp                      14.0        Exponent of the zero damping (ATM only)
    ======================== =========== ============================================

    The version of the damping can be changed after constructing the dispersion correction.
    With the `atm` boolean the three-body dispersion energy can be enabled, which is
    generally recommended.

    For a periodic cell the `ewald` dict enables the reciprocal space summation of the
    two-body dispersion energy, which removes the truncation error of the real space
    summation. An empty dict selects the default settings, the entries ``rank``,
    ``tolerance`` and ``kcut`` allow to control the accuracy. The ``mesh`` entry
    switches to the particle mesh evaluation, which is substantially faster for
    large cells. This requires a cell with
    three-dimensional periodic boundary conditions and either the rational or the zero
    damping function.

    Examples
    --------
    >>> from pyscf import gto
    >>> import dftd3.pyscf as disp
    >>> mol = gto.M(
    ...     atom='''
    ...          C   -0.189833176  -0.645396435   0.069807761
    ...          C    1.121636324  -0.354065576   0.439096514
    ...          C    1.486520953   0.962572632   0.712107225
    ...          C    0.549329390   1.989209324   0.617868956
    ...          C   -0.757627135   1.681862630   0.246856908
    ...          C   -1.138190460   0.370551816  -0.028582325
    ...          Br  -2.038462778   3.070459841   0.115165429
    ...          H    1.852935245  -1.146434699   0.514119204
    ...          H    0.825048723   3.012176989   0.829385472
    ...          H    2.502259769   1.196433556   1.000317333
    ...          H   -2.157140187   0.151608161  -0.313181471
    ...          H   -0.480820487  -1.664983631  -0.142918416
    ...          S   -4.157443472   5.729584377  -0.878761129
    ...          H   -4.823791426   4.796089466  -1.563433338
    ...          C   -2.828338520   5.970593053  -2.091189515
    ...          H   -2.167577293   6.722356639  -1.668621815
    ...          H   -2.264954814   5.054835899  -2.240198499
    ...          H   -3.218524904   6.337447714  -3.035087058
    ...          '''
    ... )
    >>> d3 = disp.DFTD3Dispersion(mol, xc="PW6B95", version="d3bj")
    >>> d3.kernel()[0]
    array(-0.01009386)
    >>> d3.version = "d3zero"  # Change to zero damping
    >>> d3.kernel()[0]
    array(-0.00574098)
    >>> d3.atm = True  # Activate three-body dispersion
    >>> d3.kernel()[0]
    array(-0.00574289)
    """

    def __init__(
        self,
        mol: Union[gto.Mole, pbc.gto.Cell],
        xc: str = "hf",
        version: str = "d3bj",
        atm: bool = False,
        param: Optional[Dict[str, float]] = None,
        ewald: Optional[Dict[str, float]] = None,
    ):
        self.mol = mol
        self.verbose = mol.verbose
        self.xc = xc
        self.param = param
        self.atm = atm
        self.version = version
        self.ewald = ewald

    def _create_model(self) -> DispersionModel:
        """Create the dispersion model for the current molecule or cell"""

        mol = self.mol

        lattice = None
        periodic = None
        if hasattr(mol, "lattice_vectors"):
            lattice = mol.lattice_vectors()
            periodic = np.array([True, True, True], dtype=bool)

        disp = DispersionModel(
            np.array([gto.charge(mol.atom_symbol(ia)) for ia in range(mol.natm)]),
            mol.atom_coords(),
            lattice=lattice,
            periodic=periodic,
        )

        if self.ewald is not None:
            # a finite system would silently fall back to real space, while
            # still using the approximated C6 coefficients
            if lattice is None:
                raise ValueError(
                    "Ewald summation requires three-dimensional periodic "
                    "boundary conditions"
                )
            disp.set_ewald_summation(
                rank=self.ewald.get("rank", 0),
                tolerance=self.ewald.get("tolerance", 0.0),
                kcut=self.ewald.get("kcut", 0.0),
                mesh=self.ewald.get("mesh", 0),
            )

        return disp

    def dump_flags(self, verbose: Optional[bool] = None) -> "DFTD3Dispersion":
        """
        Show options used for the DFT-D3 dispersion correction.
        """
        lib.logger.info(self, "** DFTD3 parameter **")
        lib.logger.info(self, "func %s", self.xc)
        lib.logger.info(
            self, "version %s", self.version + "-atm" if self.atm else self.version
        )
        return self

    def kernel(self) -> Tuple[float, np.ndarray, np.ndarray]:
        """
        Compute the DFT-D3 dispersion correction.

        The dispersion model as well as the parameters are created locally and
        not part of the state of the instance.

        Returns
        -------
        float, ndarray
            The energy and gradient of the DFT-D3 dispersion correction.

        Examples
        --------
        >>> from pyscf import gto
        >>> import dftd3.pyscf as disp
        >>> mol = gto.M(
        ...     atom='''
        ...          Br    0.000000    0.000000    1.919978
        ...          Br    0.000000    0.000000   -0.367147
        ...          N     0.000000    0.000000   -3.235006
        ...          C     0.000000    0.000000   -4.376626
        ...          H     0.000000    0.000000   -5.444276
        ...          '''
        ... )
        >>> d3 = disp.DFTD3Dispersion(mol, xc="PBE0")
        >>> energy, gradient = d3.kernel()
        >>> energy
        array(-0.00303589)
        >>> gradient
        array([[ 0.00000000e+00,  0.00000000e+00,  9.66197638e-05],
               [ 0.00000000e+00,  0.00000000e+00,  2.36000434e-04],
               [ 0.00000000e+00,  0.00000000e+00, -1.16718302e-04],
               [ 0.00000000e+00,  0.00000000e+00, -1.84332770e-04],
               [ 0.00000000e+00,  0.00000000e+00, -3.15691249e-05]])
        """
        disp = self._create_model()

        if self.param is not None:
            param = _damping_param[self.version](**self.param)
        else:
            param = _damping_param[self.version](
                method=self.xc,
                atm=self.atm,
            )

        res = disp.get_dispersion(param=param, grad=True)
        return res.get("energy"), res.get("gradient"), res.get("virial")

    def hessian(self) -> np.ndarray:
        """
        Compute the analytical second derivatives of the DFT-D3 dispersion
        correction with respect to the nuclear coordinates.

        Returns
        -------
        ndarray
            The dispersion hessian in the pyscf layout ``(natm, natm, 3, 3)``
            in Hartree per Bohr squared.

        Examples
        --------
        >>> from pyscf import gto
        >>> import dftd3.pyscf as disp
        >>> mol = gto.M(
        ...     atom='''
        ...          O  0.00000000  0.00000000 -0.73578586
        ...          H  1.44183152  0.00000000  0.36789293
        ...          H -1.44183152  0.00000000  0.36789293
        ...          '''
        ... )
        >>> d3 = disp.DFTD3Dispersion(mol, xc="PBE0")
        >>> d3.hessian().shape
        (3, 3, 3, 3)
        """
        disp = self._create_model()

        if self.param is not None:
            param = _damping_param[self.version](**self.param)
        else:
            param = _damping_param[self.version](
                method=self.xc,
                atm=self.atm,
            )

        res = disp.get_hessian(param=param)

        # (3*natm, 3*natm) with index 3*i+c -> (natm, natm, 3, 3)
        natm = self.mol.natm
        return res.get("hessian").reshape(natm, 3, natm, 3).transpose(0, 2, 1, 3)

    def reset(self, mol: Union[gto.Mole, pbc.gto.Cell]) -> "DFTD3Dispersion":
        """Reset mol and clean up relevant attributes for scanner mode"""
        self.mol = mol
        return self


class _DFTD3:
    """
    Stub class used to identify instances of the `DFTD3` class
    """

    pass


class _DFTD3Grad:
    """
    Stub class used to identify instances of the `DFTD3Grad` class
    """

    pass


class _DFTD3Hess:
    """
    Stub class used to identify instances of the `DFTD3Hess` class
    """

    pass


def d3_energy(mf: scf.hf.SCF, method: Optional[str] = None, **kwargs) -> scf.hf.SCF:
    """
    Apply DFT-D3 corrections to SCF or MCSCF methods by returning an
    instance of a new class built from the original instances class.
    The dispersion correction is stored in the `with_dftd3` attribute of
    the class.

    Parameters
    ----------
    mf: scf.hf.SCF
        The method to which DFT-D3 corrections will be applied.
    method: str, optional
        The exchange-correlation functional to use for the DFT-D3 correction.
    **kwargs
        Keyword arguments passed to the `DFTD3Dispersion` class.

    Returns
    -------
    The method with DFT-D3 corrections applied.

    Examples
    --------
    >>> from pyscf import gto, scf
    >>> import dftd3.pyscf as disp
    >>> mol = gto.M(
    ...     atom='''
    ...          N  -1.57871857  -0.04661102   0.00000000
    ...          N   1.57871857   0.04661102   0.00000000
    ...          H  -2.15862174   0.13639605   0.80956529
    ...          H  -0.84947130   0.65819321   0.00000000
    ...          H  -2.15862174   0.13639605  -0.80956529
    ...          H   2.15862174  -0.13639605  -0.80956529
    ...          H   0.84947130  -0.65819321   0.00000000
    ...          H   2.15862174  -0.13639605   0.80956529
    ...          '''
    ... )
    >>> mf = disp.d3_energy(scf.RHF(mol)).run()
    converged SCF energy = -110.932603617026
    >>> mf.kernel()
    -110.93260361702605
    """

    if not isinstance(mf, (scf.hf.SCF, mcscf.casci.CASCI)):
        raise TypeError("mf must be an instance of SCF or CASCI")

    if method is None:
        method = (
            "hf"
            if isinstance(mf, mcscf.casci.CASCI)
            else getattr(mf, "xc", "hf").lower().replace(" ", "")
        )

    with_dftd3 = DFTD3Dispersion(
        mf.mol,
        xc=method,
        **kwargs,
    )

    if isinstance(mf, _DFTD3):
        mf.with_dftd3 = with_dftd3
        return mf

    class DFTD3(_DFTD3, mf.__class__):
        def __init__(self, method: scf.hf.SCF, with_dftd3: DFTD3Dispersion):
            self.__dict__.update(method.__dict__)
            self.with_dftd3 = with_dftd3
            self._keys.update(["with_dftd3"])

        def dump_flags(self, verbose: Optional[bool] = None) -> "DFTD3":
            mf.__class__.dump_flags(self, verbose)
            if self.with_dftd3:
                self.with_dftd3.dump_flags(verbose)
            return self

        def energy_nuc(self) -> float:
            enuc = mf.__class__.energy_nuc(self)
            if self.with_dftd3:
                edisp = self.with_dftd3.kernel()[0]
                mf.scf_summary["dispersion"] = edisp
                enuc += edisp
            return enuc

        def reset(self, mol: Optional[Union[gto.Mole, pbc.gto.Cell]] = None) -> "DFTD3":
            self.with_dftd3.reset(mol)
            return mf.__class__.reset(self, mol)

        def nuc_grad_method(self) -> GradientsBase:
            scf_grad = mf.__class__.nuc_grad_method(self)
            return d3_grad(scf_grad, method=method, **kwargs)

        Gradients = lib.alias(nuc_grad_method, alias_name="Gradients")

    # only SCF and DFT classes provide a hessian method
    if hasattr(mf.__class__, "Hessian"):

        def _dftd3_hessian_method(self):
            scf_hess = mf.__class__.Hessian(self)
            return d3_hess(scf_hess, method=method, **kwargs)

        _dftd3_hessian_method.__name__ = "Hessian"
        DFTD3.Hessian = _dftd3_hessian_method

    return DFTD3(mf, with_dftd3)


def d3_grad(scf_grad: GradientsBase, **kwargs) -> GradientsBase:
    """
    Apply DFT-D3 corrections to SCF or MCSCF nuclear gradients methods
    by returning an instance of a new class built from the original class.
    The dispersion correction is stored in the `with_dftd3` attribute of
    the class.

    Parameters
    ----------
    scf_grad: rhf_grad.Gradients
        The method to which DFT-D3 corrections will be applied.
    **kwargs
        Keyword arguments passed to the `DFTD3Dispersion` class.

    Returns
    -------
    The method with DFT-D3 corrections applied.

    Examples
    --------
    >>> from pyscf import gto, scf
    >>> import dftd3.pyscf as disp
    >>> mol = gto.M(
    ...     atom='''
    ...          O  -1.65542061  -0.12330038   0.00000000
    ...          O   1.24621244   0.10268870   0.00000000
    ...          H  -0.70409026   0.03193167   0.00000000
    ...          H  -2.03867273   0.75372294   0.00000000
    ...          H   1.57598558  -0.38252146  -0.75856129
    ...          H   1.57598558  -0.38252146   0.75856129
    ...          '''
    ... )
    >>> grad = disp.d3_energy(scf.RHF(mol)).run().nuc_grad_method()
    converged SCF energy = -149.947191000075
    >>> g = grad.kernel()
    --------------- DFTD3 gradients ---------------
             x                y                z
    0 O     0.0171886976     0.0506606246     0.0000000000
    1 O     0.0383596853    -0.0459057549     0.0000000000
    2 H    -0.0313133974    -0.0125865676    -0.0000000000
    3 H     0.0066705789    -0.0380501872     0.0000000000
    4 H    -0.0154527822     0.0229409425     0.0215141991
    5 H    -0.0154527822     0.0229409425    -0.0215141991
    ----------------------------------------------
    """

    if not isinstance(scf_grad, GradientsBase):
        raise TypeError(f"scf_grad must be an instance of {GradientsBase.__name__}")

    # Ensure that the zeroth order results include DFTD3 corrections
    if not getattr(scf_grad.base, "with_dftd3", None):
        scf_grad.base = d3_energy(scf_grad.base, **kwargs)

    class DFTD3Grad(_DFTD3Grad, scf_grad.__class__):
        def grad_nuc(
            self,
            mol: Optional[Union[gto.Mole, pbc.gto.Cell]] = None,
            atmlst: Optional[np.ndarray] = None,
        ) -> np.ndarray:
            nuc_g = scf_grad.__class__.grad_nuc(self, mol, atmlst)
            with_dftd3 = getattr(self.base, "with_dftd3", None)
            if with_dftd3:
                disp_g = with_dftd3.kernel()[1]
                if atmlst is not None:
                    disp_g = disp_g[atmlst]
                nuc_g += disp_g
            return nuc_g

        def get_stress(self) -> np.ndarray:
            stress = scf_grad.__class__.get_stress(self)
            with_dftd3 = getattr(self.base, "with_dftd3", None)
            if with_dftd3:
                disp_stress = with_dftd3.kernel()[2] / abs(
                    np.linalg.det(self.cell.lattice_vectors())
                )
                stress += disp_stress
            return stress

    mfgrad = DFTD3Grad.__new__(DFTD3Grad)
    mfgrad.__dict__.update(scf_grad.__dict__)
    return mfgrad


def d3_hess(scf_hess, **kwargs):
    """
    Apply DFT-D3 corrections to SCF nuclear hessian methods by returning an
    instance of a new class built from the original class.
    The dispersion correction is stored in the `with_dftd3` attribute of
    the class.

    Parameters
    ----------
    scf_hess: hessian.rhf.Hessian
        The method to which DFT-D3 corrections will be applied.
    **kwargs
        Keyword arguments passed to the `DFTD3Dispersion` class.

    Returns
    -------
    The hessian method with DFT-D3 corrections applied.

    Examples
    --------
    >>> from pyscf import gto, scf
    >>> import dftd3.pyscf as disp
    >>> mol = gto.M(
    ...     atom='''
    ...          O  0.00000000  0.00000000 -0.73578586
    ...          H  1.44183152  0.00000000  0.36789293
    ...          H -1.44183152  0.00000000  0.36789293
    ...          '''
    ... )
    >>> mf = disp.d3_energy(scf.RHF(mol)).run()
    >>> hess = mf.Hessian().kernel()
    """

    # Ensure that the zeroth order results include DFTD3 corrections
    if not hasattr(scf_hess, "base") or not hasattr(scf_hess, "hess_nuc"):
        raise TypeError("scf_hess must be a pyscf nuclear hessian method")

    if not getattr(scf_hess.base, "with_dftd3", None):
        scf_hess.base = d3_energy(scf_hess.base, **kwargs)

    class DFTD3Hess(_DFTD3Hess, scf_hess.__class__):
        def hess_nuc(
            self,
            mol: Optional[Union[gto.Mole, pbc.gto.Cell]] = None,
            atmlst: Optional[np.ndarray] = None,
        ) -> np.ndarray:
            nuc_h = scf_hess.__class__.hess_nuc(self, mol, atmlst)
            with_dftd3 = getattr(self.base, "with_dftd3", None)
            if with_dftd3:
                disp_h = with_dftd3.hessian()
                if atmlst is not None:
                    disp_h = disp_h[np.ix_(atmlst, atmlst)]
                nuc_h += disp_h
            return nuc_h

    mfhess = DFTD3Hess.__new__(DFTD3Hess)
    mfhess.__dict__.update(scf_hess.__dict__)
    return mfhess


energy = d3_energy
"""Alias for the `d3_energy` function, which applies DFT-D3 corrections to SCF or MCSCF methods."""
grad = d3_grad
"""Alias for the `d3_grad` function, which applies DFT-D3 corrections to SCF or MCSCF nuclear gradients methods."""
hess = d3_hess
"""Alias for the `d3_hess` function, which applies DFT-D3 corrections to SCF nuclear hessian methods."""


class GeometricCounterpoiseCorrection(lib.StreamObject):
    """
    Implementation of the interface for using GCP in pyscf.

    Parameters
    ----------
    mol: gto.Mole or pbc.gto.Cell
        The molecule or periodic cell for which the GCP correction is computed.
    method: str
        The exchange-correlation functional to use for the GCP correction.
    basis: str
        The basis set to use for the GCP correction.
    """

    def __init__(
        self,
        mol: Union[gto.Mole, pbc.gto.Cell],
        method: str,
        basis: str,
    ):
        self.mol = mol
        self.verbose = mol.verbose
        self.method = method
        self.basis = basis

    def dump_flags(
        self, verbose: Optional[bool] = None
    ) -> "GeometricCounterpoiseCorrection":
        """
        Show options used for the GCP correction.
        """
        lib.logger.info(self, "** GCP parameter **")
        lib.logger.info(self, "method %s", self.method)
        lib.logger.info(self, "basis %s", self.basis)
        return self

    def kernel(self) -> Tuple[float, np.ndarray, np.ndarray]:
        mol = self.mol

        lattice = None
        periodic = None
        if hasattr(mol, "lattice_vectors"):
            lattice = mol.lattice_vectors()
            periodic = np.array([True, True, True], dtype=bool)

        gcp = GeometricCounterpoise(
            np.array([gto.charge(mol.atom_symbol(ia)) for ia in range(mol.natm)]),
            mol.atom_coords(),
            lattice=lattice,
            periodic=periodic,
            method=self.method,
            basis=self.basis,
        )

        res = gcp.get_counterpoise(grad=True)

        return res.get("energy"), res.get("gradient"), res.get("virial")

    def hessian(self) -> np.ndarray:
        """
        Compute the analytical second derivatives of the GCP correction with
        respect to the nuclear coordinates.

        Returns
        -------
        ndarray
            The counter-poise hessian in the pyscf layout ``(natm, natm, 3, 3)``
            in Hartree per Bohr squared.

        Examples
        --------
        >>> from pyscf import gto
        >>> import dftd3.pyscf as disp
        >>> mol = gto.M(
        ...     atom='''
        ...          O  0.00000000  0.00000000 -0.73578586
        ...          H  1.44183152  0.00000000  0.36789293
        ...          H -1.44183152  0.00000000  0.36789293
        ...          '''
        ... )
        >>> gcp = disp.GeometricCounterpoiseCorrection(mol, "hf3c", "minix")
        >>> gcp.hessian().shape
        (3, 3, 3, 3)
        """
        mol = self.mol

        lattice = None
        periodic = None
        if hasattr(mol, "lattice_vectors"):
            lattice = mol.lattice_vectors()
            periodic = np.array([True, True, True], dtype=bool)

        gcp = GeometricCounterpoise(
            np.array([gto.charge(mol.atom_symbol(ia)) for ia in range(mol.natm)]),
            mol.atom_coords(),
            lattice=lattice,
            periodic=periodic,
            method=self.method,
            basis=self.basis,
        )

        res = gcp.get_hessian()

        # (3*natm, 3*natm) with index 3*i+c -> (natm, natm, 3, 3)
        return (
            res.get("hessian").reshape(mol.natm, 3, mol.natm, 3).transpose(0, 2, 1, 3)
        )

    def reset(
        self, mol: Union[gto.Mole, pbc.gto.Cell]
    ) -> "GeometricCounterpoiseCorrection":
        """Reset mol and clean up relevant attributes for scanner mode"""
        self.mol = mol
        return self


class _GCP:
    pass


class _GCPGrad:
    pass


class _GCPHess:
    pass


def gcp_energy(
    mf: scf.hf.SCF, method: Optional[str] = None, basis: Optional[str] = None
) -> scf.hf.SCF:
    """
    Apply GCP corrections to SCF or MCSCF methods by returning an
    instance of a new class built from the original instances class.
    The GCP correction is stored in the `with_gcp` attribute of
    the class.

    Parameters
    ----------
    mf: scf.hf.SCF
        The method to which GCP corrections will be applied.
    method: str, optional
        The exchange-correlation functional to use for the GCP correction.
    basis: str, optional
        The basis set to use for the GCP correction.

    Returns
    -------
    The method with GCP corrections applied.
    """

    if method is None:
        method = (
            "hf"
            if isinstance(mf, mcscf.casci.CASCI)
            else getattr(mf, "xc", "hf").lower().replace(" ", "")
        )

    if basis is None:
        basis = mf.mol.basis

    with_gcp = GeometricCounterpoiseCorrection(
        mf.mol,
        method=method,
        basis=basis,
    )

    if isinstance(mf, _GCP):
        mf.with_gcp = with_gcp
        return mf

    class GCP(_GCP, mf.__class__):
        def __init__(
            self, method: scf.hf.SCF, with_gcp: GeometricCounterpoiseCorrection
        ):
            self.__dict__.update(method.__dict__)
            self.with_gcp = with_gcp
            self._keys.update(["with_gcp"])

        def dump_flags(self, verbose: Optional[bool] = None) -> "GCP":
            mf.__class__.dump_flags(self, verbose)
            if self.with_gcp:
                self.with_gcp.dump_flags(verbose)
            return self

        def energy_nuc(self) -> float:
            enuc = mf.__class__.energy_nuc(self)
            if self.with_gcp:
                edisp = self.with_gcp.kernel()[0]
                mf.scf_summary["dispersion"] = edisp
                enuc += edisp
            return enuc

        def reset(self, mol: Optional[Union[gto.Mole, pbc.gto.Cell]] = None) -> "GCP":
            self.with_gcp.reset(mol)
            return mf.__class__.reset(self, mol)

        def nuc_grad_method(self) -> GradientsBase:
            scf_grad = mf.__class__.nuc_grad_method(self)
            return gcp_grad(scf_grad, method=method, basis=basis)

        Gradients = lib.alias(nuc_grad_method, alias_name="Gradients")

    # only SCF and DFT classes provide a hessian method
    if hasattr(mf.__class__, "Hessian"):

        def _gcp_hessian_method(self):
            scf_hess = mf.__class__.Hessian(self)
            return gcp_hess(scf_hess, method=method, basis=basis)

        _gcp_hessian_method.__name__ = "Hessian"
        GCP.Hessian = _gcp_hessian_method

    return GCP(mf, with_gcp)


def gcp_grad(scf_grad: GradientsBase, **kwargs) -> GradientsBase:
    """
    Apply GCP corrections to SCF or MCSCF nuclear gradients methods
    by returning an instance of a new class built from the original class.
    The GCP correction is stored in the `with_gcp` attribute of
    the class.

    Parameters
    ----------
    scf_grad: rhf_grad.Gradients
        The method to which GCP corrections will be applied.
    **kwargs
        Keyword arguments passed to the `GeometricCounterpoiseCorrection` class.

    Returns
    -------
    The method with GCP corrections applied.
    """
    # Ensure that the zeroth order results include GCP corrections
    if not getattr(scf_grad.base, "with_gcp", None):
        scf_grad.base = gcp_energy(scf_grad.base, **kwargs)

    class GCPGrad(_GCPGrad, scf_grad.__class__):
        def grad_nuc(
            self,
            mol: Optional[Union[gto.Mole, pbc.gto.Cell]] = None,
            atmlst: Optional[np.ndarray] = None,
        ) -> np.ndarray:
            nuc_g = scf_grad.__class__.grad_nuc(self, mol, atmlst)
            with_gcp = getattr(self.base, "with_gcp", None)
            if with_gcp:
                gcp_g = with_gcp.kernel()[1]
                if atmlst is not None:
                    gcp_g = gcp_g[atmlst]
                nuc_g += gcp_g
            return nuc_g

        def get_stress(self) -> np.ndarray:
            stress = scf_grad.__class__.get_stress(self)
            with_gcp = getattr(self.base, "with_gcp", None)
            if with_gcp:
                disp_stress = with_gcp.kernel()[2] / abs(
                    np.linalg.det(self.cell.lattice_vectors())
                )
                stress += disp_stress
            return stress

    mfgrad = GCPGrad.__new__(GCPGrad)
    mfgrad.__dict__.update(scf_grad.__dict__)
    return mfgrad


def gcp_hess(scf_hess, **kwargs):
    """
    Apply GCP corrections to SCF nuclear hessian methods by returning an
    instance of a new class built from the original class.
    The GCP correction is stored in the `with_gcp` attribute of
    the class.

    Parameters
    ----------
    scf_hess: hessian.rhf.Hessian
        The method to which GCP corrections will be applied.
    **kwargs
        Keyword arguments passed to the `GeometricCounterpoiseCorrection` class.

    Returns
    -------
    The hessian method with GCP corrections applied.

    Examples
    --------
    >>> from pyscf import gto, scf
    >>> import dftd3.pyscf as disp
    >>> mol = gto.M(
    ...     atom='''
    ...          O  0.00000000  0.00000000 -0.73578586
    ...          H  1.44183152  0.00000000  0.36789293
    ...          H -1.44183152  0.00000000  0.36789293
    ...          '''
    ... )
    >>> mf = disp.gcp_energy(scf.RHF(mol)).run()
    >>> hess = mf.Hessian().kernel()
    """

    if not hasattr(scf_hess, "base") or not hasattr(scf_hess, "hess_nuc"):
        raise TypeError("scf_hess must be a pyscf nuclear hessian method")

    # Ensure that the zeroth order results include GCP corrections
    if not getattr(scf_hess.base, "with_gcp", None):
        scf_hess.base = gcp_energy(scf_hess.base, **kwargs)

    class GCPHess(_GCPHess, scf_hess.__class__):
        def hess_nuc(
            self,
            mol: Optional[Union[gto.Mole, pbc.gto.Cell]] = None,
            atmlst: Optional[np.ndarray] = None,
        ) -> np.ndarray:
            nuc_h = scf_hess.__class__.hess_nuc(self, mol, atmlst)
            with_gcp = getattr(self.base, "with_gcp", None)
            if with_gcp:
                gcp_h = with_gcp.hessian()
                if atmlst is not None:
                    gcp_h = gcp_h[np.ix_(atmlst, atmlst)]
                nuc_h += gcp_h
            return nuc_h

    mfhess = GCPHess.__new__(GCPHess)
    mfhess.__dict__.update(scf_hess.__dict__)
    return mfhess
