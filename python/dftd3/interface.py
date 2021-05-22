# This file is part of s-dftd3.
# SPDX-Identifier: LGPL-3.0-or-later
#
# s-dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# s-dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.
"""Wrapper around the C-API of the s-dftd3 shared library."""

from typing import List, Optional, Union
import numpy as np


from .libdftd3 import (
    ffi as _ffi,
    lib as _lib,
    new_error,
    handle_error,
    new_structure,
    new_d3_model,
    new_zero_damping,
    load_zero_damping,
    new_rational_damping,
    load_rational_damping,
    new_mzero_damping,
    load_mzero_damping,
    new_mrational_damping,
    load_mrational_damping,
)


class Structure:
    """
    .. Molecular structure data

    Represents a wrapped structure object in ``s-dftd3``.
    The molecular structure data object has a fixed number of atoms
    and immutable atomic identifiers
    """

    _mol = _ffi.NULL

    def __init__(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
    ):
        """Create new molecular structure data"""
        if positions.size % 3 != 0:
            raise ValueError("Expected tripels of cartesian coordinates")

        if 3 * numbers.size != positions.size:
            raise ValueError("Dimension missmatch between numbers and positions")

        self._natoms = len(numbers)
        _numbers = np.ascontiguousarray(numbers, dtype="i4")
        _positions = np.ascontiguousarray(positions, dtype=float)

        if lattice is not None:
            if lattice.size != 9:
                raise ValueError("Invalid lattice provided")
            _lattice = np.ascontiguousarray(lattice, dtype="float")
        else:
            _lattice = None

        if periodic is not None:
            if periodic.size != 3:
                raise ValueError("Invalid periodicity provided")
            _periodic = np.ascontiguousarray(periodic, dtype="bool")
        else:
            _periodic = None

        self._mol = new_structure(
            self._natoms,
            _cast("int*", _numbers),
            _cast("double*", _positions),
            _cast("double*", _lattice),
            _cast("bool*", _periodic),
        )

    def __len__(self):
        return self._natoms

    def update(
        self,
        positions: np.ndarray,
        lattice: Optional[np.ndarray] = None,
    ) -> None:
        """Update coordinates and lattice parameters, both provided in
        atomic units (Bohr).
        The lattice update is optional also for periodic structures.

        Generally, only the cartesian coordinates and the lattice parameters
        can be updated, every other modification, regarding total charge,
        total spin, boundary condition, atomic types or number of atoms
        requires the complete reconstruction of the object.
        """

        if 3 * len(self) != positions.size:
            raise ValueError("Dimension missmatch for positions")
        _positions = np.ascontiguousarray(positions, dtype="float")

        if lattice is not None:
            if lattice.size != 9:
                raise ValueError("Invalid lattice provided")
            _lattice = np.ascontiguousarray(lattice, dtype="float")
        else:
            _lattice = None

        _error = new_error()
        _lib.dftd3_update_structure(
            _error,
            self._mol,
            _cast("double*", _positions),
            _cast("double*", _lattice),
        )

        handle_error(_error)


class DampingParam:
    """Abstract base class for damping parameters"""

    _param = _ffi.NULL

    def __init__(self, **kwargs):
        """Create new damping parameter from method name or explicit data"""

        if "method" in kwargs and kwargs["method"] is not None:
            self._param = self.load_param(**kwargs)
        else:
            self._param = self.new_param(**kwargs)

    @staticmethod
    def load_param(method, **kwargs):
        raise NotImplementedError("Child class has to define parameter loading")

    @staticmethod
    def new_param(**kwargs):
        raise NotImplementedError("Child class has to define parameter construction")


class RationalDampingParam(DampingParam):
    def __init__(self, **kwargs):
        DampingParam.__init__(self, **kwargs)

    @staticmethod
    def load_param(method, **kwargs):
        _method = _ffi.new("char[]", method.encode())
        return load_rational_damping(
            _method,
            kwargs.get("s9", 1.0) > 0.0,
        )

    @staticmethod
    def new_param(**kwargs):
        try:
            return new_rational_damping(
                kwargs.get("s6", 1.0),
                kwargs["s8"],
                kwargs.get("s9", 1.0),
                kwargs["a1"],
                kwargs["a2"],
                kwargs.get("alp", 14.0),
            )
        except KeyError as e:
            raise RuntimeError("Constructor requires argument for " + str(e))


class ZeroDampingParam(DampingParam):
    def __init__(self, **kwargs):
        DampingParam.__init__(self, **kwargs)

    @staticmethod
    def load_param(method, **kwargs):
        _method = _ffi.new("char[]", method.encode())
        return load_zero_damping(
            _method,
            kwargs.get("s9", 1.0) > 0.0,
        )

    @staticmethod
    def new_param(**kwargs):
        try:
            return new_zero_damping(
                kwargs.get("s6", 1.0),
                kwargs["s8"],
                kwargs.get("s9", 1.0),
                kwargs["rs6"],
                kwargs.get("rs8", 1.0),
                kwargs.get("alp", 14.0),
            )
        except KeyError as e:
            raise RuntimeError("Constructor requires argument for " + str(e))


class ModifiedRationalDampingParam(DampingParam):
    def __init__(self, **kwargs):
        DampingParam.__init__(self, **kwargs)

    @staticmethod
    def load_param(method, **kwargs):
        _method = _ffi.new("char[]", method.encode())
        return load_mrational_damping(
            _method,
            kwargs.get("s9", 1.0) > 0.0,
        )

    @staticmethod
    def new_param(**kwargs):
        try:
            return new_mrational_damping(
                kwargs.get("s6", 1.0),
                kwargs["s8"],
                kwargs.get("s9", 1.0),
                kwargs["a1"],
                kwargs["a2"],
                kwargs.get("alp", 14.0),
            )
        except KeyError as e:
            raise RuntimeError("Constructor requires argument for " + str(e))


class ModifiedZeroDampingParam(DampingParam):
    def __init__(self, **kwargs):
        DampingParam.__init__(self, **kwargs)

    @staticmethod
    def load_param(method, **kwargs):
        _method = _ffi.new("char[]", method.encode())
        return load_mzero_damping(
            _method,
            kwargs.get("s9", 1.0) > 0.0,
        )

    @staticmethod
    def new_param(**kwargs):
        try:
            return new_mzero_damping(
                kwargs.get("s6", 1.0),
                kwargs["s8"],
                kwargs.get("s9", 1.0),
                kwargs["rs6"],
                kwargs.get("rs8", 1.0),
                kwargs.get("alp", 14.0),
                kwargs["bet"],
            )
        except KeyError as e:
            raise RuntimeError("Constructor requires argument for " + str(e))


class DispersionModel(Structure):
    """
    .. Dispersion model
    """

    _disp = _ffi.NULL

    def __init__(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
    ):
        Structure.__init__(self, numbers, positions, lattice, periodic)

        self._disp = new_d3_model(self._mol)

    def get_dispersion(self, param: DampingParam, grad: bool) -> dict:
        """Perform actual evaluation of the dispersion correction"""

        _error = new_error()
        _energy = _ffi.new("double *")
        if grad:
            _gradient = np.zeros((len(self), 3))
            _sigma = np.zeros((3, 3))
        else:
            _gradient = None
            _sigma = None

        _lib.dftd3_get_dispersion(
            _error,
            self._mol,
            self._disp,
            param._param,
            _energy,
            _cast("double*", _gradient),
            _cast("double*", _sigma),
        )

        handle_error(_error)

        results = dict(energy=_energy[0])
        if _gradient is not None:
            results.update(gradient=_gradient)
        if _sigma is not None:
            results.update(virial=_sigma)
        return results


def _cast(ctype, array):
    """Cast a numpy array to a FFI pointer"""
    return _ffi.NULL if array is None else _ffi.cast(ctype, array.ctypes.data)


def _ref(ctype, value):
    """Create a reference to a value"""
    if value is None:
        return _ffi.NULL
    ref = _ffi.new(ctype + "*")
    ref[0] = value
    return ref
