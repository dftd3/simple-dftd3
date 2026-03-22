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
QCSchema Support
================

Integration with the `QCArchive infrastructure <http://docs.qcarchive.molssi.org>`_.

If the QCElemental package is installed the ``dftd3.qcschema`` module becomes
importable and provides the ``run_qcschema`` function supporting QCSchema v1.
If the QCElemental package is >=0.50.0, ``dftd3.qcschema`` supports QCSchema v1
and v2, returning whichever version was submitted. Note that Python 3.14+ only
works with QCSchema v2 due to Pydantic restrictions.

This module provides a way to translate QCSchema or QCElemental Atomic Input
into a format understandable by the ``dftd3`` API which in turn provides the
calculation results in a QCSchema compatible format.

Supported keywords are

======================== =========== ============================================
 Keyword                  Default     Description
======================== =========== ============================================
 level_hint               None        Dispersion correction level
 params_tweaks            None        Optional dict with the damping parameters
 pair_resolved            False       Enable pairwise resolved dispersion energy
======================== =========== ============================================

Allowed level hints are ``"d3bj"``, ``"d3zero"``, ``"d3bjm"``/``"d3mbj"``,
``"d3mzero"``/``"d3zerom"``, ``"d3op"``, and ``"d3cso"``.

The params_tweaks dict contains the damping parameters, at least s8, a1 and a2
must be provided for rational damping, while s8 and rs6 are required in case
of zero damping. For CSO damping, a1 must be provided.

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

.. note::

    input_data.model.method with a full method name and input_data.keywords["params_tweaks"]
    cannot be provided at the same time. It is an error to provide both options at the
    same time.

Example
-------

>>> from dftd3.qcschema import run_qcschema
>>> import qcelemental as qcel
>>> atomic_input = qcel.models.AtomicInput(
...     molecule = qcel.models.Molecule(
...         symbols = ["O", "H", "H"],
...         geometry = [
...             0.00000000000000,  0.00000000000000, -0.73578586109551,
...             1.44183152868459,  0.00000000000000,  0.36789293054775,
...            -1.44183152868459,  0.00000000000000,  0.36789293054775
...         ],
...     ),
...     driver = "energy",
...     model = {
...         "method": "TPSS-D3(BJ)",
...     },
...     keywords = {},
... )
...
>>> atomic_result = run_qcschema(atomic_input)
>>> atomic_result.return_result
-0.00042042440936212056
"""

from typing import Union
from .interface import (
    DispersionModel,
    RationalDampingParam,
    ZeroDampingParam,
    ModifiedRationalDampingParam,
    ModifiedZeroDampingParam,
    OptimizedPowerDampingParam,
    CSODampingParam,
)
from .library import get_api_version
import numpy as np
import qcelemental as qcel


_supported_drivers = [
    "energy",
    "gradient",
]

_available_levels = [
    "d3bj",
    "d3zero",
    "d3bjm",
    "d3mbj",
    "d3zerom",
    "d3mzero",
    "d3op",
    "d3cso",
]

_damping_param = {
    "d3bj": RationalDampingParam,
    "d3zero": ZeroDampingParam,
    "d3bjm": ModifiedRationalDampingParam,
    "d3mbj": ModifiedRationalDampingParam,
    "d3zerom": ModifiedZeroDampingParam,
    "d3mzero": ModifiedZeroDampingParam,
    "d3op": OptimizedPowerDampingParam,
    "d3cso": CSODampingParam,
}

_clean_dashlevel = str.maketrans("", "", "()")


def run_qcschema(
    input_data: Union[dict, qcel.models.AtomicInput, "qcel.models.v2.AtomicInput"]
) -> Union[qcel.models.AtomicResult, "qcel.models.v2.AtomicResult"]:
    """Perform disperson correction based on an atomic inputmodel"""

    v2_available = hasattr(qcel.models, "v2")

    if v2_available and isinstance(input_data, qcel.models.v2.AtomicInput):
        atomic_input = input_data
    elif isinstance(input_data, qcel.models.AtomicInput):
        atomic_input = input_data
    elif v2_available and input_data.get("specification"):
        atomic_input = qcel.models.v2.AtomicInput(**input_data)
    else:
        atomic_input = qcel.models.AtomicInput(**input_data)

    if (schver := atomic_input.schema_version) == 1:
        from qcelemental.models import AtomicInput, AtomicResult, ComputeError

        ret_data = atomic_input.dict()
    elif schver == 2:
        from qcelemental.models.v2 import AtomicInput, AtomicResult, ComputeError, FailedOperation

        ret_data = {"input_data": atomic_input, "extras": {}, "molecule": atomic_input.molecule}

    provenance = {
        "creator": "s-dftd3",
        "version": get_api_version(),
        "routine": "dftd3.qcschema.run_qcschema",
    }
    success = False
    return_result = 0.0
    properties = {}

    # Since it is a level hint we a forgiving if it is not present,
    # we are much less forgiving if the wrong level is hinted here.
    atin_keywords = atomic_input.keywords if schver == 1 else atomic_input.specification.keywords
    _level = atin_keywords.get("level_hint", "d3bj")
    if _level.lower() not in _available_levels:
        error=ComputeError(
            error_type="input error",
            error_message="Level '{}' is invalid for this dispersion correction".format(
                _level
            ),
        )
        if schver == 1:
            ret_data.update(
                provenance=provenance,
                success=success,
                properties=properties,
                return_result=return_result,
                error=error,
            )
            return AtomicResult(**ret_data)
        elif schver == 2:
            return FailedOperation(input_data=atomic_input, error=error)

    # Check if the method is provided and strip the “dashlevel” from the method
    _method = atomic_input.model.method.split("-") if schver == 1 else atomic_input.specification.model.method.split("-")
    if _method[-1].lower().translate(_clean_dashlevel) == _level.lower():
        _method.pop()
    _method = "-".join(_method)
    if len(_method) == 0:
        _method = None

    # Obtain the parameters for the damping function
    _input_param = atin_keywords.get("params_tweaks", {"method": _method})

    try:
        param = _damping_param[_level](
            **_input_param,
        )

        disp = DispersionModel(
            atomic_input.molecule.atomic_numbers[atomic_input.molecule.real],
            atomic_input.molecule.geometry[atomic_input.molecule.real],
        )

        driver = atomic_input.driver if schver == 1 else atomic_input.specification.driver
        res = disp.get_dispersion(
            param=param,
            grad=driver == "gradient",
        )
        extras = {"dftd3": res}

        if driver == "gradient":
            if all(atomic_input.molecule.real):
                fullgrad = res.get("gradient")
            else:
                ireal = np.argwhere(atomic_input.molecule.real).reshape((-1))
                fullgrad = np.zeros_like(atomic_input.molecule.geometry)
                fullgrad[ireal, :] = res.get("gradient")

        properties.update(return_energy=res.get("energy"))

        if atin_keywords.get("pair_resolved", False):
            res = disp.get_pairwise_dispersion(param=param)
            extras["dftd3"].update(res)

        success = driver in _supported_drivers
        if driver == "energy":
            return_result = properties["return_energy"]
        elif driver == "gradient":
            return_result = fullgrad
        else:
            ret_data.update(
                error=ComputeError(
                    error_type="input error",
                    error_message="Calculation succeeded but invalid driver request provided",
                ),
            )

        ret_data["extras"].update(extras)

    except (RuntimeError, TypeError) as e:
        ret_data.update(
            error=ComputeError(
                error_type="input error", error_message=str(e)
            ),
        ),

    ret_data.update(
        provenance=provenance,
        success=success,
        properties=properties,
        return_result=return_result,
    )

    if schver == 2 and "error" in ret_data:
        return FailedOperation(input_data=atomic_input, error=ret_data["error"])

    return AtomicResult(**ret_data)
