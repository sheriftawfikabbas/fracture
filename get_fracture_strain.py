import glob
import ntpath
import numpy as np
import pandas as pd

results_file = open("results_revised.csv", "w")

folders = glob.glob("gitlog_cnt_fracture/*")

# folders = glob.glob("sample_structures/test_code/*")


def read(f):
    ff = open(f)
    t = ff.readlines()
    ff.close()
    t = t[2:]
    t = [[float(k) for k in x.split()[1:]] for x in t]
    return np.array(t)


def differences(a):
    a = sorted(a)
    return np.array(a[1 : len(a) - 1]) - np.array(a[0 : len(a) - 2])


def check(dm_initial, dm_current):
    for i in range(len(dm_current)):
        if dm_current[i] > 4:
            print("fracture!")
            # print a snapshot of the two structures, current and previous
            return True
    return False


for folder in folders:

    results_file.write(folder + ",")
    structures = glob.glob(folder + "/fracture_*.xyz")
    folder_bottom = folder.split("/")[1]
    print(folder_bottom)

    ids = sorted(
        [int(ntpath.basename(x).split("_")[2].replace(".xyz", "")) for x in structures]
    )

    prefix = ntpath.basename(structures[0]).split("_")
    prefix = prefix[0] + "_" + prefix[1] + "_"
    a_initial = read(folder + "/" + prefix + str(0) + ".xyz")
    dm_initial = differences(a_initial[:, 2])
    found = False
    height_initial = max(a_initial[:, 2]) - min(a_initial[:, 2])
    for im in ids[1:]:
        a_current = read(folder + "/" + prefix + str(im) + ".xyz")
        height_current = max(a_current[:, 2]) - min(a_current[:, 2])
        s = (height_current - height_initial) / height_initial * 100
        dm_current = differences(a_current[:, 2])
        # print(folder, im)
        if check(dm_initial, dm_current):
            results_file.write(folder + "," + str(im) + "," + str(s) + ",\n")
            results_file.flush()
            found = True
            break
    if not found:
        results_file.write(folder + "," + str(im) + "," + str(s) + ",nofracture\n")
        results_file.flush()
results_file.close()
