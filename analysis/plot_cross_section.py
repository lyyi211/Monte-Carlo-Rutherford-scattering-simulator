import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent

data_file = (
    project_dir
    / "data"
    / "monte_carlo_scattering50000.csv"
)

results_dir = project_dir / "results"
results_dir.mkdir(exist_ok=True)


# Parameters

K = 227.514312
E = 5.0         
b_max = 100.0   


# Read num scattering angles

theta_deg = []

with open(data_file, "r") as file:
    reader = csv.DictReader(file)

    for row in reader:
        theta_deg.append(
            float(row["theta_numerical"])
        )

theta_deg = np.array(theta_deg)

N = len(theta_deg)

print("Number of Monte Carlo events:", N)


# Angular bins

# Using 5-degree bins. Simulation does not contain angles below about 25.6 deg
# because b is limited to b_max = 100 fm.

bin_width = 5.0

theta_min_rad = 2.0 * np.arctan(
    K / (2.0 * E * b_max)
)

theta_min_deg = np.degrees(theta_min_rad)

bin_edges_deg = np.concatenate(
    (
        [theta_min_deg],
        np.arange(30.0, 181.0, bin_width)
    )
)


# Count events in each angular bin
counts, edges_deg = np.histogram(
    theta_deg,
    bins=bin_edges_deg
)


# Bin centers
theta_centers_deg = (
    edges_deg[:-1] +
    edges_deg[1:]
) / 2.0


# Convert edges and centers to radians
edges_rad = np.radians(edges_deg)

theta_centers_rad = np.radians(
    theta_centers_deg
)


# Cross section represented by 1 event
sigma_total_sampled = (
    np.pi * b_max**2
)

sigma_per_event = (
    sigma_total_sampled / N
)


# Cross section in each bin
delta_sigma = (
    counts * sigma_per_event
)


# Solid angle of each bin
delta_omega = (
    2.0
    * np.pi
    * (
        np.cos(edges_rad[:-1])
        -
        np.cos(edges_rad[1:])
    )
)


# Differential cross section
dsigma_domega_mc = (
    delta_sigma /
    delta_omega
)


# Statistical uncertainty
dsigma_error = (
    np.sqrt(counts)
    * sigma_per_event
    / delta_omega
)


# Rutherford analytical cross section
theta_curve_deg = np.linspace(
    theta_min_deg,
    180.0,
    1000
)

theta_curve_rad = np.radians(
    theta_curve_deg
)


rutherford = (
    (K / (4.0 * E))**2
    /
    np.sin(theta_curve_rad / 2.0)**4
)


# Remove empty bins

mask = counts > 0

theta_plot = theta_centers_deg[mask]

mc_plot = dsigma_domega_mc[mask]

error_plot = dsigma_error[mask]


# PLOT 1

plt.figure(figsize=(8, 5))

plt.errorbar(
    theta_plot,
    mc_plot,
    yerr=error_plot,
    fmt="o",
    capsize=3,
    label="Monte Carlo + RK4"
)

plt.plot(
    theta_curve_deg,
    rutherford,
    label="Rutherford theory"
)


plt.yscale("log")


plt.xlabel(
    "Scattering angle $\\theta$ (degrees)"
)

plt.ylabel(
    r"$d\sigma/d\Omega$ (log scale)(fm$^2$/sr)"
)

plt.title(
    "Rutherford differential scattering cross section"
)

plt.grid(
    True,
    which="both",
    alpha=0.3
)

plt.legend()

plt.tight_layout()


plt.savefig(
    results_dir
    / "rutherford_differential_cross_section.png",
    dpi=300
)

plt.show()


#=============================================================================
# Rutherford theory witn bins

A = (K / (4.0 * E))**2

theta1 = edges_rad[:-1]
theta2 = edges_rad[1:]


# Cross section integrated over each angular bin

delta_sigma_rutherford = (
    4.0
    * np.pi
    * A
    * (
        1.0 / np.sin(theta1 / 2.0)**2
        -
        1.0 / np.sin(theta2 / 2.0)**2
    )
)


rutherford_bin_average = (
    delta_sigma_rutherford
    / delta_omega
)


rutherford_bin_average_plot = (
    rutherford_bin_average[mask]
)


# Error

relative_residual = (
    mc_plot
    -
    rutherford_bin_average_plot
) / rutherford_bin_average_plot


relative_residual_error = (
    error_plot
    /
    rutherford_bin_average_plot
)


relative_residual_percent = (
    100.0 * relative_residual
)

relative_residual_error_percent = (
    100.0 * relative_residual_error
)



print("\nRelative residuals:\n")

for theta, residual, uncertainty in zip(
    theta_plot,
    relative_residual_percent,
    relative_residual_error_percent
):
    print(
        f"theta = {theta:6.2f} deg   "
        f"residual = {residual:8.3f} %   "
        f"stat. uncertainty = {uncertainty:8.3f} %"
    )

# Chi-squared goodness of fit

chi_squared = np.sum(
    (relative_residual / relative_residual_error)**2
)

degrees_of_freedom = len(relative_residual)

reduced_chi_squared = (
    chi_squared / degrees_of_freedom
)

print("\nGoodness of fit:")
print(f"chi^2 = {chi_squared:.3f}")
print(f"degrees of freedom = {degrees_of_freedom}")
print(f"reduced chi^2 = {reduced_chi_squared:.3f}")


# PLOT 2

plt.figure(figsize=(8, 5))

plt.errorbar(
    theta_plot,
    relative_residual_percent,
    yerr=relative_residual_error_percent,
    fmt="o",
    capsize=3,
    label="Monte Carlo + RK4"
)


plt.axhline(
    y=0.0,
    linestyle="--",
    label="Rutherford theory"
)


plt.xlabel(
    "Scattering angle $\\theta$ (degrees)"
)

plt.ylabel(
    "Relative residual (%)"
)

plt.title(
    "Relative deviation from Rutherford theory"
)

plt.grid(
    True,
    alpha=0.3
)

plt.legend()

plt.tight_layout()


plt.savefig(
    results_dir
    / "rutherford_relative_residuals.png",
    dpi=300
)

plt.show()
