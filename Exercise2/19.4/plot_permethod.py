import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

methods = {
    "RK2 (explicit)"    : "ex19_4_rk2.csv",
    "RK4 (explicit)"    : "ex19_4_rk4.csv",
    "IRK Gauss-2"       : "ex19_4_irk_gauss2.csv",
}

for label, filename in methods.items():
    df = pd.read_csv(filename)
    taus = sorted(df["tau"].unique())

    plt.figure(figsize=(10,5))
    for tau in taus:
        sub = df[df["tau"] == tau]
        plt.plot(sub["t"], sub["y"], label=f"y(t), tau={tau}")
        plt.plot(sub["t"], sub["v"], linestyle="--", label=f"v(t), tau={tau}")
    plt.xlabel("t")
    plt.ylabel("y(t) and v(t)")
    plt.title(f"Time evolution ({label})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"time_{label.replace(' ','_')}.png", dpi=200)

    plt.figure(figsize=(6,6))
    for tau in taus:
        sub = df[df["tau"] == tau]
        plt.plot(sub["y"], sub["v"], label=f"tau={tau}")
    plt.xlabel("y")
    plt.ylabel("v")
    plt.title(f"Phase plot ({label})")
    plt.grid(True)
    plt.axis("equal")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"phase_{label.replace(' ','_')}.png", dpi=200)

plt.show()
