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

from typing import Iterator

import numpy as np
import pytest
from pytest import approx, raises

try:
    import ase
    from dftd3.ase import DFTD3
    from ase.build import molecule
    from ase.calculators.emt import EMT
except ModuleNotFoundError:
    ase = None


def get_calcs(calc) -> Iterator[ase.calculators.calculator.Calculator]:
    if hasattr(calc, "mixer"):
        calc = calc.mixer
    yield from calc.calcs


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_ase_scand4():
    thr = 1.0e-6

    forces = np.array(
        [
            [-0.00000000e-00, -0.00000000e-00, -6.83426991e-05],
            [-0.00000000e-00, +3.44839555e-04, +7.21176947e-04],
            [+6.80565391e-22, -3.44839555e-04, +7.21176947e-04],
            [+2.12676685e-23, +3.26671388e-20, -1.60555514e-03],
            [+3.21599856e-04, +5.54947267e-04, +6.86106874e-04],
            [-3.21599856e-04, +5.54947267e-04, +6.86106874e-04],
            [-3.21599856e-04, -5.54947267e-04, +6.86106874e-04],
            [+3.21599856e-04, -5.54947267e-04, +6.86106874e-04],
            [+1.87155483e-21, +2.87678390e-04, -1.25644177e-03],
            [-3.40282696e-22, -2.87678390e-04, -1.25644177e-03],
        ]
    )

    atoms = molecule("methylenecyclopropane")
    atoms.calc = DFTD3(method="SCAN", damping="d3bj")

    assert atoms.get_potential_energy() == approx(-0.03880921894019244, abs=thr)
    assert atoms.get_forces() == approx(forces, abs=thr)

    atoms.calc = DFTD3(method="SCAN", damping="d3bj").add_calculator(EMT())
    assert atoms.get_potential_energy() == approx(3.6452960962398406, abs=thr)
    energies = [calc.get_potential_energy() for calc in get_calcs(atoms.calc)]
    assert energies == approx([-0.03880921894019244, 3.684105315180033], abs=thr)


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_ase_tpssd4():
    thr = 1.0e-6

    forces = np.array(
        [
            [+1.21727790e-03, +1.98579200e-03, -1.16371697e-02],
            [-5.82484114e-04, +9.01770290e-03, +7.78537640e-03],
            [-4.30031958e-03, +4.63213536e-03, -4.56657109e-03],
            [-1.16941383e-03, -8.39071556e-03, +1.60593512e-02],
            [-6.90354443e-03, -5.07801933e-03, -1.75396161e-03],
            [+1.03561818e-02, -1.68908740e-02, -2.74225314e-03],
            [+5.59001294e-03, +3.35129491e-03, -9.24429928e-04],
            [-5.13316989e-03, +6.07626858e-03, +3.89454026e-05],
            [+3.35952011e-03, +3.95424504e-03, -5.65438002e-04],
            [-2.13140242e-03, +2.77295425e-03, -4.76829804e-04],
            [+4.33961724e-03, -1.51731003e-03, -7.01598391e-04],
            [-4.64227572e-03, +8.65258554e-05, -5.15421318e-04],
        ]
    )

    atoms = molecule("C2H6CHOH")
    atoms.calc = DFTD3(method="TPSS", damping="d3zero")

    assert atoms.get_potential_energy() == approx(-0.14230914516094673, abs=thr)
    assert atoms.get_forces() == approx(forces, abs=thr)

    atoms.calc = DFTD3(method="TPSS", damping="d3zero").add_calculator(EMT())
    assert atoms.get_potential_energy() == approx(4.963774668847532, abs=thr)
    energies = [calc.get_potential_energy() for calc in get_calcs(atoms.calc)]
    assert energies == approx([-0.14230914516094673, 5.106083814008478], abs=thr)
