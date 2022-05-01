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

import numpy as np
import pytest
from pytest import approx, raises

try:
    import pyscf
    from pyscf import lib, gto, scf
    import dftd3.pyscf as disp
except ModuleNotFoundError:
    pyscf = None


@pytest.mark.skipif(pyscf is None, reason="requires pyscf")
def test_energy_r2scan_d3():

    mol = gto.M(
        atom="""
             C   -0.755422531  -0.796459123  -1.023590391
             C    0.634274834  -0.880017014  -1.075233285
             C    1.406955202   0.199695367  -0.653144334
             C    0.798863737   1.361204515  -0.180597909
             C   -0.593166787   1.434312023  -0.133597923
             C   -1.376239198   0.359205222  -0.553258516
             I   -1.514344238   3.173268101   0.573601106
             H    1.110906949  -1.778801728  -1.440619836
             H    1.399172302   2.197767355   0.147412751
             H    2.486417780   0.142466525  -0.689380574
             H   -2.454252250   0.422581120  -0.512807958
             H   -1.362353593  -1.630564523  -1.348743149
             S   -3.112683203   6.289227834   1.226984439
             H   -4.328789697   5.797771251   0.973373089
             C   -2.689135032   6.703163830  -0.489062886
             H   -1.684433029   7.115457372  -0.460265708
             H   -2.683867206   5.816530502  -1.115183775
             H   -3.365330613   7.451201412  -0.890098894
             """
    )

    d3 = disp.DFTD3Dispersion(mol, xc="r2SCAN")
    assert d3.kernel()[0] == approx(-0.00578401192369041, abs=1.0e-7)


@pytest.mark.skipif(pyscf is None, reason="requires pyscf")
def test_gradient_b97m_d3():
    mol = gto.M(
        atom="""
             H    0.002144194   0.361043475   0.029799709
             C    0.015020592   0.274789738   1.107648016
             C    1.227632658   0.296655040   1.794629427
             C    1.243958826   0.183702791   3.183703934
             C    0.047958213   0.048915002   3.886484583
             C   -1.165135654   0.026954348   3.200213281
             C   -1.181832083   0.139828643   1.810376587
             H    2.155807907   0.399177037   1.249441585
             H    2.184979344   0.198598553   3.716170761
             H    0.060934662  -0.040672756   4.964014252
             H   -2.093220602  -0.078628959   3.745125056
             H   -2.122845437   0.123257119   1.277645797
             Br  -0.268325907  -3.194209024   1.994458950
             C    0.049999933  -5.089197474   1.929391171
             F    0.078949601  -5.512441335   0.671851563
             F    1.211983937  -5.383996300   2.498664481
             F   -0.909987405  -5.743747328   2.570721738
             """
    )
    ref = np.array(
        [
            [+7.13721248e-07, +2.19571763e-05, -3.77372946e-05],
            [+9.19838860e-07, +3.53459763e-05, -1.43306994e-06],
            [+7.43860881e-06, +3.78237447e-05, +8.46031238e-07],
            [+8.06120927e-06, +3.79834948e-05, +8.58427570e-06],
            [+1.16592466e-06, +3.62585085e-05, +1.16326308e-05],
            [-3.69381337e-06, +3.39047971e-05, +6.92483428e-06],
            [-3.05404225e-06, +3.29484247e-05, +1.80766271e-06],
            [+3.51228183e-05, +2.08136972e-05, -1.76546837e-05],
            [+3.49762054e-05, +1.66544908e-05, +2.14435772e-05],
            [+1.57516340e-06, +1.41373959e-05, +4.21574793e-05],
            [-3.35392428e-05, +1.49030766e-05, +2.29976305e-05],
            [-3.38817253e-05, +1.82002569e-05, -1.72487448e-05],
            [-2.15610724e-05, -1.87935101e-04, -3.02815495e-05],
            [+1.27580963e-06, -5.96841724e-05, -5.99713166e-06],
            [+9.01173808e-07, -2.23010304e-05, -7.96228701e-06],
            [+7.42062176e-06, -2.79631452e-05, +7.03703317e-07],
            [-3.84119900e-06, -2.30475903e-05, +1.21693625e-06],
        ]
    )

    d3 = disp.DFTD3Dispersion(mol, xc="r2SCAN")
    assert d3.kernel()[1] == approx(ref, abs=1.0e-7)


@pytest.mark.skipif(pyscf is None, reason="requires pyscf")
def test_energy_hf():
    mol = gto.M(
        atom="""
             N  -1.57871857  -0.04661102   0.00000000
             N   1.57871857   0.04661102   0.00000000
             H  -2.15862174   0.13639605   0.80956529
             H  -0.84947130   0.65819321   0.00000000
             H  -2.15862174   0.13639605  -0.80956529
             H   2.15862174  -0.13639605  -0.80956529
             H   0.84947130  -0.65819321   0.00000000
             H   2.15862174  -0.13639605   0.80956529
             """
    )
    mf = disp.energy(scf.RHF(mol))
    assert mf.kernel() == approx(-110.93260361702605, abs=1.0e-8)


@pytest.mark.skipif(pyscf is None, reason="requires pyscf")
def test_gradient_hf():
    mol = gto.M(
        atom="""
             O  -1.65542061  -0.12330038   0.00000000
             O   1.24621244   0.10268870   0.00000000
             H  -0.70409026   0.03193167   0.00000000
             H  -2.03867273   0.75372294   0.00000000
             H   1.57598558  -0.38252146  -0.75856129
             H   1.57598558  -0.38252146   0.75856129
             """
    )
    ref = np.array(
        [
            [+1.71886976e-02, +5.06606246e-02, +3.88042826e-16],
            [+3.83596853e-02, -4.59057549e-02, -1.90819582e-16],
            [-3.13133974e-02, -1.25865676e-02, -2.79182059e-16],
            [+6.67057892e-03, -3.80501872e-02, -9.99808574e-17],
            [-1.54527822e-02, +2.29409425e-02, +2.15141991e-02],
            [-1.54527822e-02, +2.29409425e-02, -2.15141991e-02],
        ]
    )
    grad = disp.energy(scf.RHF(mol)).run().nuc_grad_method()
    assert grad.kernel() == approx(ref, abs=1.0e-8)
