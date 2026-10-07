import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent

data_file = (
    project_dir
    / "data"
    / "monte_carlo_impact_parameters.csv"
)

results_dir = project_dir / "results"
results_dir.mkdir(exist_ok=True)


# Read data

b_values = []

with open(data_file, "r") as file:
    reader = csv.DictReader(file)

    for row in reader:
        b_values.append(float(row["b"]))
        

# Parameters

b_max = 100.0

b_theory = np.linspace(0.0, b_max, 500)

p_theory = 2.0 * b_theory / (b_max ** 2)


# Histogram 

plt.figure(figsize=(8, 5))

plt.hist(
    b_values,
    bins=50,
    density=True,
    alpha=0.6,
    label="Monte Carlo simulation"
)

plt.plot(
    b_theory,
    p_theory,
    linewidth=2,
    label="Theoretical distribution"
)

plt.xlabel("Impact parameter b (fm)")
plt.ylabel("Probability density (1/fm)")

plt.title(
    "Monte Carlo sampling of the impact parameter"
)

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(
    results_dir / "impact_parameter_distribution.png",
    dpi=300
)

plt.show()
