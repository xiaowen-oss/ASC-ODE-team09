import numpy as np
import matplotlib.pyplot as plt
import csv
import os

# Save location (inside your repo)
save_path = os.path.join(os.getcwd(), "Exercise2", "phase_ex18_5.png")

theta = []
omega = []

with open("Exercise2/ex18_5.csv", "r") as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        theta.append(float(row[1]))
        omega.append(float(row[2]))

theta = np.array(theta)
omega = np.array(omega)

plt.figure(figsize=(7,5))
plt.plot(theta, omega, linewidth=1.5)

plt.xlabel("θ (theta)")
plt.ylabel("ω (omega)")
plt.title("Phase Plot of the Pendulum")
plt.grid(True)
plt.tight_layout()

plt.savefig(save_path, dpi=200)
plt.show()

print("Saved to:", save_path)
