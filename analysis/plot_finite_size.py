
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


project_dir = Path(__file__).resolve().parent.parent

data_file = project_dir / "data" / "finite_size_scattering.csv"
results_dir = project_dir / "results"

results_dir.mkdir(exist_ok=True)

data = pd.read_csv(data_file)


# Constants

eSquared = 1.439964
Z_ALPHA = 2.0
Z_GOLD = 79.0

K = eSquared * Z_ALPHA * Z_GOLD

A_GOLD = 197.0
r0 = 1.2

R = r0 * A_GOLD**(1.0 / 3.0)

U_surface = K / R
U_center = 3.0 * K / (2.0 * R)


# PLOT 1
# Scattering angle

energies = [5, 35, 40, 50, 75, 100]

fig, axes = plt.subplots(2, 3, figsize=(14, 8))

for ax, E in zip(axes.flat, energies):

    subset = data[data["E"] == E].sort_values("b")

    ax.plot(
        subset["b"],
        subset["theta_point"],
        "-",
        linewidth=2,
        label="Point nucleus"
    )

    ax.plot(
        subset["b"],
        subset["theta_sphere"],
        "--",
        linewidth=2,
        label="Uniform sphere"
    )

    ax.set_title(f"E = {E} MeV")
    ax.set_xlabel("Impact parameter b (fm)")
    ax.set_ylabel("Scattering angle (degrees)")

    ax.set_xlim(0, 15)
    ax.set_ylim(0, 185)

    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

fig.suptitle(
    "Scattering angle: point nucleus vs uniform sphere",
    fontsize=14
)

plt.tight_layout()

plt.savefig(
    results_dir / "finite_size_scattering_angles.png",
    dpi=300
)

plt.show()


# PLOT 2
# Finite size deviation heatmap

pivot = data.pivot(
    index="E",
    columns="b",
    values="delta_theta"
)

energies_grid = pivot.index.to_numpy(dtype=float)
b_grid = pivot.columns.to_numpy(dtype=float)

fig, ax = plt.subplots(figsize=(10, 6))

im = ax.pcolormesh(
    b_grid,
    energies_grid,
    pivot.values,
    shading="nearest",
    cmap="coolwarm",
    vmin=-180,
    vmax=0
)

ax.set_xlim(0, 15)
ax.set_ylim(5, 100)

ax.set_xlabel("Impact parameter b (fm)")
ax.set_ylabel("Alpha particle energy E (MeV)")

ax.set_title(
    "Finite-size deviation from point-nucleus scattering"
)

cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Delta theta (degrees)")

plt.tight_layout()

plt.savefig(
    results_dir / "finite_size_deviation_heatmap.png",
    dpi=300
)

plt.show()


# PLOT 3
# Potential energy

r = np.linspace(0.1, 30.0, 1000)

U_point = K / r

U_sphere = np.where(
    r < R,
    K / (2.0 * R) * (3.0 - r**2 / R**2),
    K / r
)

plt.figure(figsize=(9, 6))

plt.plot(
    r,
    U_point,
    label="Point nucleus"
)

plt.plot(
    r,
    U_sphere,
    label="Uniform sphere"
)

plt.axvline(
    R,
    linestyle="--",
    label=f"Nuclear radius R = {R:.2f} fm"
)

plt.axhline(
    U_surface,
    linestyle=":",
    label=f"U(R) = {U_surface:.2f} MeV"
)

plt.axhline(
    U_center,
    linestyle=":",
    label=f"U(0) = {U_center:.2f} MeV"
)

plt.xlim(0, 30)
plt.ylim(0, 120)

plt.xlabel("Distance r (fm)")
plt.ylabel("Potential energy U(r) (MeV)")

plt.title(
    "Coulomb potential energy: point vs finite-size nucleus"
)

plt.grid(True, alpha=0.3)
plt.legend()

plt.tight_layout()

plt.savefig(
    results_dir / "finite_size_potential_energy.png",
    dpi=300
)

plt.show()
