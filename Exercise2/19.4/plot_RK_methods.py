import numpy as np
import matplotlib.pyplot as plt

def load_csv(filename, tau_value):
    data = np.genfromtxt(filename, delimiter=",", names=True)
    mask = np.isclose(data["tau"], tau_value)
    return data["t"][mask], data["y"][mask], data["v"][mask]

methods = {
    "RK2 (explicit)"  : "ex19_4_rk2.csv",
    "RK4 (explicit)"  : "ex19_4_rk4.csv",
    "IRK Gauss-2"     : "ex19_4_irk_gauss2.csv",
}

taus = [0.1, 0.05, 0.01]

for tau in taus:
    plt.figure(figsize=(10,5))
    for label, filename in methods.items():
        t, y, v = load_csv(filename, tau)
        plt.plot(t, y, label=label)
    plt.xlabel("t")
    plt.ylabel("y")
    plt.title(f"Time evolution (tau={tau})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"time_evolution_tau_{tau}.png", dpi=200)

    plt.figure(figsize=(6,6))
    for label, filename in methods.items():
        t, y, v = load_csv(filename, tau)
        plt.plot(y, v, label=label)
    plt.xlabel("y")
    plt.ylabel("v")
    plt.title(f"Phase plot (tau={tau})")
    plt.grid(True)
    plt.axis("equal")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"phase_plot_tau_{tau}.png", dpi=200)

plt.show()
