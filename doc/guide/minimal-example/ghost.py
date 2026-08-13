import numpy as np

from dftd3.interface import DispersionModel, RationalDampingParam

num = np.array([8, 1, 1, 8, 1, 1])
xyz = np.array(
    [
        [-2.5142928288858868, 1.8879764891638844, 0.0000000000000000],
        [-2.0098096791094657, 3.6303285583904779, 0.0000000000000000],
        [-0.9545622395768059, 0.92963588958819932, 0.0000000000000000],
        [1.8995798491583649, -1.4196911493711153, 0.0000000000000000],
        [1.7895424492078418, -2.5141248888855077, -1.4467274293585810],
        [1.7895424492078418, -2.5141248888855077, 1.4467274293585810],
    ]
)

model = DispersionModel(num, xyz)
model.set_ghost_atoms([3, 4, 5])  # 0-based atom indices of the second fragment
res = model.get_dispersion(RationalDampingParam(method="PBE0"), grad=False)
print(f"Dispersion energy with ghost atoms is {res['energy']:13.10f} Hartree")
