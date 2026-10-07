import csv
from pathlib import Path
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent.parent

data_file = ROOT / "data" / "scattering_angles.csv"

results_dir = ROOT / "results"
results_dir.mkdir(exist_ok=True)


# Read data

b = []
theta_numerical = []
theta_theory = []
relative_error = []

with open(data_file, "r") as file:
    reader = csv.DictReader(file)

    for row in reader:
        b.append(float(row["b"]))
        theta_numerical.append(float(row["theta_numerical"]))
        theta_theory.append(float(row["theta_theory"]))
        relative_error.append(
            float(row["relative_error"]) * 100.0
        )


# PLOT 1
# Numerical RK4 vs analytical Rutherford solution

plt.figure(figsize=(8, 5))

plt.plot(
    b,
    theta_theory,
    label="Rutherford analytical solution",
    linewidth=2
)

plt.scatter(
    b,
    theta_numerical,
    label="Numerical RK4 simulation",
    s=35
)

plt.xlabel("Impact parameter b (fm)")
plt.ylabel("Scattering angle θ (degrees)")

plt.title(
    "Rutherford scattering: numerical vs analytical solution"
)

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(
    ROOT / "results" / "scattering_angle_comparison.png",
    dpi=300
)

plt.show()


# PLOT 2
# Relative numerical error

plt.figure(figsize=(8, 5))

plt.plot(
    b,
    relative_error,
    marker="o",
    linewidth=1.5
)

plt.xlabel("Impact parameter b (fm)")
plt.ylabel("Relative error (%)")

plt.title(
    "Relative error of numerical scattering angle"
)

plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig(
    ROOT / "results" / "relative_error_num_angle.png",
    dpi=300
)

plt.show()
