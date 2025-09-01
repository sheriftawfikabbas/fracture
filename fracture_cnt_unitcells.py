from ase.build import nanotube
from oganesson import OgStructure
import os
import sys
import torch
import glob, ntpath

torch.set_default_device("cuda")
model = sys.argv[1]
step = float(sys.argv[2])
method = sys.argv[3]

total_c = 200

g = glob.glob("cnt_unitcells/*.cif")
for gg in g:
    try:
        print(">> Working on", gg)
        cif_file = ntpath.basename(gg).replace(".cif", "")
        folder = (
            "cnt_fracture/"
            + model
            + "_"
            + str(step)
            + "_"
            + method
            + "_gpu_"
            + cif_file
        )
        if not os.path.isdir(folder):
            os.mkdir(folder)
        og = OgStructure(file_name=gg)
        c = og().lattice.c
        og.make_supercell([1, 1, int(total_c / c)])
        og().to(folder + "/" + cif_file + "_unoptimised.cif")
        og.relax(relax_cell=False)
        og().to(folder + "/" + cif_file + "_optimised.cif")
        og.fracture(
            1.25,
            write_intermediate=True,
            write_intermediate_file_extension="xyz",
            method=method,
            intermediates_folder=folder + "/",
            model=model,
            translation_step=step,
            fmax=0.1,
            relaxation_steps=500,
        )
        og.to_ase().write(folder + "/" + cif_file + ".xyz")
    except:
        continue
