import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read CSV
df = pd.read_csv("crank_nicolson.csv")

taus = sorted(df["tau"].unique())

# Time evolution plot

plt.figure(figsize=(12,5))
for tau in taus:
    data = df[df["tau"] == tau]
    t = data["t"]
    y = data["y"]
    v = data["v"]

    plt.plot(t, y, label=f"x(t), tau={tau}")
    plt.plot(t, v, '--', label=f"v(t), tau={tau}")

plt.xlabel("t")
plt.ylabel("x(t) and v(t)")
plt.title("Time Evolution (Crank–Nicolson)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# Phase plot

plt.figure(figsize=(6,6))
for tau in taus:
    data = df[df["tau"] == tau]
    y = data["y"]
    v = data["v"]

    plt.plot(y, v, label=f"tau={tau}")

plt.xlabel("x")
plt.ylabel("v")
plt.title("Phase Plot (Crank–Nicolson)")
plt.grid(True)
plt.axis("equal")
plt.legend()
plt.tight_layout()
plt.show()
