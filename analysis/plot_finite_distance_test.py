import csv
from pathlib import Path
import matplotlib.pyplot as plt


script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent

data_file = project_dir / "data" / "finite_distance_test.csv"

results_dir = project_dir / "results"
results_dir.mkdir(parents=True, exist_ok=True)


# Read data

x0 = []
theta_numerical = []
theta_theory = []
absolute_error = []

with open(data_file, "r") as file:
    reader = csv.DictReader(file)

    for row in reader:
        x0.append(float(row["x0"]))
        theta_numerical.append(float(row["theta_numerical"]))
        theta_theory.append(float(row["theta_theory"]))
        absolute_error.append(
            abs(float(row["absolute_error"]))
        )


# PLOT 1
# Finite-distance convergence

plt.figure(figsize=(8, 5))

plt.plot(
    x0,
    theta_numerical,
    marker="o",
    label="Numerical RK4"
)

plt.axhline(
    y=theta_theory[0],
    linestyle="--",
    label="Rutherford analytical solution"
)

plt.xlabel("Starting distance $x_0$ (fm)")
plt.ylabel("Scattering angle $\\theta$ (degrees)")

plt.title("Convergence with starting distance")

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(
    results_dir / "finite_distance_convergence.png",
    dpi=300
)

plt.show()


# PLOT 2
# Finite-distance error

plt.figure(figsize=(8, 5))

plt.loglog(
    x0,
    absolute_error,
    marker="o",
    label="Numerical error"
)



reference = [
    absolute_error[0] * x0[0] / x
    for x in x0
]

plt.loglog(
    x0,
    reference,
    linestyle="--",
    label=r"$\propto 1/x_0$"
)

plt.xlabel("Starting distance $x_0$ (fm)")
plt.ylabel(
    r"$|\theta_{\rm num}-\theta_{\rm theory}|$ (degrees)"
)

plt.title("Finite-distance error (log-log scale)")


plt.grid(
    True,
    which="major",
    alpha=0.3
)

plt.legend()
plt.tight_layout()

plt.savefig(
    results_dir / "finite_distance_error.png",
    dpi=300
)

plt.show()
