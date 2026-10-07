import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent

data_file = project_dir / "data" / "monte_carlo_scattering.csv"

results_dir = project_dir / "results"
results_dir.mkdir(exist_ok=True)


# Read data

theta_numerical = []
theta_theory = []

with open(data_file, "r") as file:
    reader = csv.DictReader(file)

    for row in reader:
        theta_numerical.append(
            float(row["theta_numerical"])
        )

        theta_theory.append(
            float(row["theta_theory"])
        )

theta_numerical = np.array(theta_numerical)
theta_theory = np.array(theta_theory)


# Histogram 1
# Scattering angle distribution 

plt.figure(figsize=(8, 5))

plt.hist(
    theta_numerical,
    bins=50,
    alpha=0.7
)

plt.xlabel("Scattering angle θ (degrees)")
plt.ylabel("Number of events")

plt.title(
    "Monte Carlo distribution of scattering angles"
)

plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig(
    results_dir / "monte_carlo_scattering_angle_distribution.png",
    dpi=300
)

plt.show()


# Histogram 2
# Angle error

angle_error = theta_numerical - theta_theory

plt.figure(figsize=(8, 5))

plt.hist(
    angle_error,
    bins=50,
    alpha=0.7
)

plt.xlabel("θ numerical − θ theory (degrees)")
plt.ylabel("Number of events")

plt.title(
    "Numerical error distribution"
)

plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig(
    results_dir / "monte_carlo_angle_error_distribution.png",
    dpi=300
)

plt.show()
