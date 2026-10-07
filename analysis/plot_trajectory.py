from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent

data = pd.read_csv(
    ROOT / "data" / "single_particle_trajectory_b30.csv"
)

x = data["x"]
y = data["y"]

# distance from the nucleus
r = np.sqrt(x**2 + y**2)

# point of closest approach
i_min = r.idxmin()

x_closest = x[i_min]
y_closest = y[i_min]
r_min = r[i_min]


plt.figure(figsize=(8, 6))

plt.plot(x, y, label="Alpha particle trajectory for b=30fm")
plt.scatter(0, 0, marker="o", label="Gold nucleus")
plt.scatter(x_closest, y_closest, marker="x",
            label="Closest approach")

plt.xlabel("x [fm]")
plt.ylabel("y [fm]")
plt.title("Rutherford single-particle scattering near the nucleus")

plt.grid(True)
plt.legend()


# zoom near the scattering region
plt.xlim(-200, 150)
plt.ylim(-50, 200)


plt.tight_layout()
plt.savefig(
    ROOT / "results" / "single_particle_trajectory_b30.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()
plt.show()
