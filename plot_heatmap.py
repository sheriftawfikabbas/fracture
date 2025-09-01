import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm

plt.rcParams.update({"font.size": 25})


window_min = 14
window_max = 19

# Implementation of Eq. 1 in the paper

def get_theta(m, n):
    return abs(np.arctan(-np.sqrt(3) * n / (n + 2 * m)))


def get_d(a, m, n):
    return np.sqrt(3) * a / np.pi * np.sqrt(n**2 + n * m + m**2)


def f(theta, v):
    return 1 / 2 * ((1 - v) + (1 + v) * np.cos(2 * theta))


def strenth(a, m, n):
    C = 55
    alpha = 0.5
    v = 0.16
    return C / f(get_theta(m, n), v) * get_d(a, m, n) ** (-alpha)


def cmap2rgb(cmap, step):
    v = getattr(cm, cmap)(step, bytes=True)
    return [x / 255 for x in v]


B = pd.read_csv("fracture_strain.csv")
x = B.m.unique()
y = B.n.unique()

x.sort()
y.sort()

s = B.s

s_table = []
s_dict = {}
for yi in y:
    rowX = []
    for xi in x:
        if (xi == 2 and yi == 0) or xi < yi:
            continue
        rowX += [B.loc[(B.m == xi) & (B.n == yi)].s.values[0]]
        s_dict[(xi, yi)] = B.loc[(B.m == xi) & (B.n == yi)].s.values[0]
    s_table += [rowX]


d_axis = [i for i in range(2, 19)]
theta_axis = [i for i in range(0, 19)]
d_axis.reverse()

B_table_df = pd.DataFrame(s_table, columns=d_axis, index=theta_axis)

plt.figure(figsize=(10, 10))
ax = sns.heatmap(B_table_df, xticklabels=d_axis, yticklabels=theta_axis, square=True)
plt.xlabel("$m$")
plt.ylabel("$n$")

plt.savefig("results", bbox_inches="tight", transparent=False, dpi=300)
plt.clf()

# Normalise
B.s -= B.s.min()


d_axis = np.linspace(0, 40, 500)
theta_axis = np.linspace(0, 40, 500)
z = np.ndarray([500, 500])
my_z = np.ndarray([500, 500])
marks_x = []
marks_y = []

experimental_mn = [
    (13, 12),
    (14, 12),
    (14, 11),
    (15, 8),
    (15, 12),
    (16, 14),
    (17, 9),
]

experimental_dtheta = []
for i, r in B.iterrows():
    s = strenth(1.41, r.m, r.n)
    theta = get_theta(r.m, r.n) * 180 / np.pi
    d = get_d(1.41, r.m, r.n)

    ix_min = 0
    d_v = abs(d - d_axis[ix_min])
    for ix in range(len(d_axis)):
        x = d_axis[ix]
        if abs(d - x) < d_v:
            d_v = abs(d - x)
            ix_min = ix
    iy_min = 0
    theta_v = abs(theta - theta_axis[iy_min])
    for iy in range(len(theta_axis)):
        y = theta_axis[iy]
        if abs(theta - y) < theta_v:
            theta_v = abs(theta - y)
            iy_min = iy
    z[ix_min, iy_min] = s
    my_z[ix_min, iy_min] = r.s
    marks_x += [ix_min]
    marks_y += [iy_min]
    if (r.m, r.n) in experimental_mn:
        experimental_dtheta += [(ix_min, iy_min)]


d_every = []
d_axis_every = []

for i in range(0, len(d_axis), 50):
    d_axis_every += [i]
    d_every += [round(d_axis[i])]

for i in range(0, len(d_axis), 50):
    if d_axis[i] > 10:
        xlim_min = i
        break
for i in range(0, len(d_axis), 50):
    if d_axis[i] > 22:
        xlim_max = i
        break

theta_every = []
theta_axis_every = []
for i in range(0, len(theta_axis), 50):
    theta_axis_every += [i]
    theta_every += [round(theta_axis[i])]


for i in range(0, len(theta_axis), 50):
    if theta_axis[i] > 15:
        ylim_min = i
        break
for i in range(0, len(theta_axis), 50):
    if theta_axis[i] > 32:
        ylim_max = i
        break

fig = plt.figure(figsize=(10, 10))
for i in range(500):
    for j in range(500):
        if z[i][j] != 0:
            if (i, j) in experimental_dtheta:
                circle = plt.Circle((i, j), 10, color=(0, 0, 0), lw=1)
                fig.gca().add_artist(circle)
            circle = plt.Circle(
                (i, j),
                4,
                color=cmap2rgb(
                    "Spectral", (1 - (z[i][j] - window_min) / (window_max - window_min))
                ),
            )
            fig.gca().add_artist(circle)

plt.xlabel("Diameter $(\AA)$")
plt.ylabel("Chiral angle ($^{o}$)")
plt.xticks(d_axis_every, labels=d_every)
plt.yticks(theta_axis_every, labels=theta_every)
plt.savefig("results_d_theta", bbox_inches="tight", transparent=False, dpi=300)

plt.clf()


fig = plt.figure(figsize=(10, 10))
for i in range(500):
    for j in range(500):
        if my_z[i][j] != 0:
            if (i, j) in experimental_dtheta:
                circle = plt.Circle((i, j), 10, color=(0, 0, 0), lw=1)
                fig.gca().add_artist(circle)
            circle = plt.Circle(
                (i, j),
                4,
                color=cmap2rgb(
                    "Spectral",
                    (
                        1
                        - (my_z[i][j] / my_z.max() * (window_max - window_min))
                        / (window_max - window_min)
                    ),
                ),
            )
            fig.gca().add_artist(circle)

print(z.min(), my_z.min())
print(z.max(), my_z.max())
plt.xlabel("Diameter $(\AA)$")
plt.ylabel("Chiral angle ($^{o}$)")
plt.xticks(d_axis_every, labels=d_every)
plt.yticks(theta_axis_every, labels=theta_every)
plt.savefig("results_my_d_theta", bbox_inches="tight", transparent=False, dpi=300)

plt.clf()

fig, axes = plt.subplots()
fig.set_figheight(10)
fig.set_figwidth(10)
for i in range(499, -1, -1):
    for j in range(499, -1, -1):
        if z[i][j] != 0:
            if (i, j) in experimental_dtheta:
                circle = plt.Circle(
                    (i, j),
                    radius=(z[i][j] - window_min) / (window_max - window_min) * 30,
                    color=cmap2rgb(
                        "Spectral",
                        (1 - (z[i][j] - window_min) / (window_max - window_min)),
                    ),
                )
                fig.gca().add_artist(circle)
                circle = plt.Circle(
                    (i, j),
                    radius=(z[i][j] - window_min) / (window_max - window_min) * 30,
                    fill=False,
                )
                axes.set_aspect(1)
                fig.gca().add_artist(circle)

plt.xlabel("Diameter $(\AA)$")
plt.ylabel("Chiral angle ($^{o}$)")
plt.xticks(d_axis_every, labels=d_every)
plt.yticks(theta_axis_every, labels=theta_every)
plt.xlim([xlim_min, xlim_max])
plt.ylim([ylim_min, ylim_max])
plt.savefig("results_d_theta_window", bbox_inches="tight", transparent=False, dpi=300)

plt.clf()


R = 30
fig, axes = plt.subplots()
fig.set_figheight(10)
fig.set_figwidth(10)
for i in range(500):
    for j in range(500):
        if my_z[i][j] != 0:
            if (i, j) in experimental_dtheta:
                circle = plt.Circle(
                    (i, j),
                    radius=(my_z[i][j] / my_z.max() * (window_max - window_min))
                    / (window_max - window_min)
                    * R,
                    color=cmap2rgb(
                        "Spectral",
                        (
                            1
                            - (my_z[i][j] / my_z.max() * (window_max - window_min))
                            / (window_max - window_min)
                        ),
                    ),
                )
                fig.gca().add_artist(circle)
                circle = plt.Circle(
                    (i, j),
                    radius=(my_z[i][j] / my_z.max() * (window_max - window_min))
                    / (window_max - window_min)
                    * R,
                    fill=False,
                )
                axes.set_aspect(1)
                fig.gca().add_artist(circle)

plt.xlabel("Diameter $(\AA)$")
plt.ylabel("Chiral angle ($^{o}$)")
plt.xticks(d_axis_every, labels=d_every)
plt.yticks(theta_axis_every, labels=theta_every)
plt.xlim([xlim_min, xlim_max])
plt.ylim([ylim_min, ylim_max])
plt.savefig(
    "results_my_d_theta_window", bbox_inches="tight", transparent=False, dpi=300
)

plt.clf()
