from ase.build import nanotube
from oganesson import OgStructure
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import torch

torch.set_default_device("cuda")


def y_fitted(x):
    return (
        coeff[0]
        + x * coeff[1]
        + x**2 * coeff[2]
        + x**3 * coeff[3]
        + x**4 * coeff[4]
        + x**5 * coeff[5]
        + x**6 * coeff[6]
    )


for i in range(0, 18 + 1):
    for j in range(0, i + 1):
        energies = []
        for k in range(-10, 10 + 1):
            atoms = nanotube(i, j, length=1, bond=1.42, symbol="C", vacuum=20)
            atoms.set_pbc = True
            og = OgStructure(atoms)
            og = og.scale([1, 1, 1 + k / 100])
            print(og().lattice)
            og.relax(relax_cell=False)
            energies += [[k, og.total_energy]]
        f = np.array(energies)
        x = f[:, 0]
        y = f[:, 1]

        # Use polyfit to get the coefficient of a 6-degree polynomial
        coeff = np.polynomial.polynomial.polyfit(x, y, 6)

        # Then find its roots using polyroots
        roots = np.polynomial.polynomial.polyroots(
            (
                coeff[1],
                2 * coeff[2],
                3 * coeff[3],
                4 * coeff[4],
                5 * coeff[5],
                6 * coeff[6],
            )
        )
        print(roots)

        y_fitted_vector = y_fitted(x)

        print(y_fitted(roots))

        # Then create a comparisong plot
        plt.figure(figsize=(10, 10))
        plt.plot(x, y)
        plt.plot(x, y_fitted_vector)
        plt.ylabel("y")
        plt.xlabel("x")
        plt.savefig(
            "cnt_unitcells_fits/cnt_" + str(i) + "_" + str(j), bbox_inches="tight"
        )

        for x in range(len(roots)):
            if x <= 10 and x >= -10:
                atoms = nanotube(i, j, length=1, bond=1.42, symbol="C", vacuum=20)
                atoms.set_pbc = True
                og = OgStructure(atoms)
                og = og.scale([1, 1, 1 + x / 100])

                og().to("cnt_unitcells/cnt_" + str(i) + "_" + str(j) + "_opt.cif")
