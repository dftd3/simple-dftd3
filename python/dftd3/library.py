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
"""
Thin wrapper around the CFFI extension. This module mainly acts as a guard
for importing the libdftd3 extension and also provides some FFI based wrappers
for memory handling.
"""

import functools

try:
    from ._libdftd3 import ffi, lib
except ImportError:
    raise ImportError("DFT-D3 C extension unimportable, cannot use C-API")


def get_api_version() -> str:
    """Return the current API version from s-dftd3.
    For easy usage in C the API version is provided as

    10000 * major + 100 * minor + patch

    For Python we want something that looks like a semantic version again.
    """
    api_version = lib.dftd3_get_version()
    return "{}.{}.{}".format(
        api_version // 10000,
        api_version % 10000 // 100,
        api_version % 100,
    )


def _delete_error(mol):
    """Delete a DFT-D3 error handle object"""
    ptr = ffi.new("dftd3_error *")
    ptr[0] = mol
    lib.dftd3_delete_error(ptr)


def new_error():
    """Create new DFT-D3 error handle object"""
    return ffi.gc(lib.dftd3_new_error(), _delete_error)


def error_check(func):
    """Handle errors for library functions that require an error handle"""

    @functools.wraps(func)
    def handle_error(*args, **kwargs):
        """Run function and than compare context"""
        _err = new_error()
        value = func(_err, *args, **kwargs)
        if lib.dftd3_check_error(_err):
            _message = ffi.new("char[]", 512)
            lib.dftd3_get_error(_err, _message, ffi.NULL)
            raise RuntimeError(ffi.string(_message).decode())
        return value

    return handle_error


def _delete_structure(mol):
    """Delete a DFT-D3 molecular structure data object"""
    ptr = ffi.new("dftd3_structure *")
    ptr[0] = mol
    lib.dftd3_delete_structure(ptr)


def new_structure(natoms, numbers, positions, lattice, periodic):
    """Create new molecular structure data"""
    return ffi.gc(
        error_check(lib.dftd3_new_structure)(
            natoms,
            numbers,
            positions,
            lattice,
            periodic,
        ),
        _delete_structure,
    )


def _delete_model(disp):
    """Delete a DFT-D3 dispersion model object"""
    ptr = ffi.new("dftd3_model *")
    ptr[0] = disp
    lib.dftd3_delete_model(ptr)


def new_d3_model(mol):
    """Create new D3 dispersion model"""
    return ffi.gc(error_check(lib.dftd3_new_d3_model)(mol), _delete_model)


set_model_realspace_cutoff = error_check(lib.dftd3_set_model_realspace_cutoff)


def _delete_param(param):
    """Delete a DFT-D3 damping parameteter object"""
    ptr = ffi.new("dftd3_param *")
    ptr[0] = param
    lib.dftd3_delete_param(ptr)


def new_zero_damping(s6, s8, s9, rs6, rs8, alp):
    """Create new zero damping parameters"""
    return ffi.gc(
        error_check(lib.dftd3_new_zero_damping)(s6, s8, s9, rs6, rs8, alp),
        _delete_param,
    )


def load_zero_damping(method, atm):
    """Load zero damping parameters from internal storage"""
    return ffi.gc(error_check(lib.dftd3_load_zero_damping)(method, atm), _delete_param)


def new_rational_damping(s6, s8, s9, a1, a2, alp):
    """Create new rational damping parameters"""
    return ffi.gc(
        error_check(lib.dftd3_new_rational_damping)(s6, s8, s9, a1, a2, alp),
        _delete_param,
    )


def load_rational_damping(method, atm):
    """Load rational damping parameters from internal storage"""
    return ffi.gc(
        error_check(lib.dftd3_load_rational_damping)(method, atm), _delete_param
    )


def new_mzero_damping(s6, s8, s9, rs6, rs8, alp, bet):
    """Create new modified zero damping parameters"""
    return ffi.gc(
        error_check(lib.dftd3_new_mzero_damping)(s6, s8, s9, rs6, rs8, alp, bet),
        _delete_param,
    )


def load_mzero_damping(method, atm):
    """Load modified zero damping parameters from internal storage"""
    return ffi.gc(error_check(lib.dftd3_load_mzero_damping)(method, atm), _delete_param)


def new_mrational_damping(s6, s8, s9, a1, a2, alp):
    """Create new modified rational damping parameters"""
    return ffi.gc(
        error_check(lib.dftd3_new_mrational_damping)(s6, s8, s9, a1, a2, alp),
        _delete_param,
    )


def load_mrational_damping(method, atm):
    """Load modified rational damping parameters from internal storage"""
    return ffi.gc(
        error_check(lib.dftd3_load_mrational_damping)(method, atm), _delete_param
    )


def new_optimizedpower_damping(s6, s8, s9, a1, a2, alp, bet):
    """Create new optimized power damping parameters"""
    return ffi.gc(
        error_check(lib.dftd3_new_optimizedpower_damping)(s6, s8, s9, a1, a2, alp, bet),
        _delete_param,
    )


def load_optimizedpower_damping(method, atm):
    """Load optimized power damping parameters from internal storage"""
    return ffi.gc(
        error_check(lib.dftd3_load_optimizedpower_damping)(method, atm), _delete_param
    )


update_structure = error_check(lib.dftd3_update_structure)
get_dispersion = error_check(lib.dftd3_get_dispersion)
get_pairwise_dispersion = error_check(lib.dftd3_get_pairwise_dispersion)


def _ref(ctype, value):
    """Create a reference to a value"""
    if value is None:
        return ffi.NULL
    ref = ffi.new(ctype + "*")
    ref[0] = value
    return ref
