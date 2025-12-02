import numpy as np
import matplotlib.pyplot as plt

def load_csv(filename, tau_value):
    data = np.genfromtxt(filename, delimiter=",", names=True)
    mask = np.isclose(data["tau"], tau_value)
    return data["t"][mask], data["y"][mask], data["v"][mask]

methods = {
    "RK2 (explicit)"        : "ex19_4_rk2.csv",
    "RK4 (explicit)"        : "ex19_4_rk4.csv",
    "IRK Gauss-2"           : "ex19_4_irk_gauss2.csv",
    "Radau IIA (implicit)"  : "ex19_4_radau2a.csv",
}


taus = [0.01, 0.1, 0.5, 1.0, 2.0]

for tau in taus:
    plt.figure(figsize=(8,4))
    for label, filename in methods.items():
        t, y, v = load_csv(filename, tau)

        # exact solution for harmonic oscillator
        y_exact = np.cos(t)
        v_exact = -np.sin(t)

        # error in phase space
        err = np.sqrt((y - y_exact)**2 + (v - v_exact)**2)

        plt.plot(t, err, label=label)

    plt.xlabel("t")
    plt.ylabel("error norm")
    plt.title(f"Error vs time (tau={tau})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"error_tau_{tau}.png", dpi=200)

plt.show()
