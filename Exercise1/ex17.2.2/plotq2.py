import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#input data
df = pd.read_csv("improved_euler.csv")

# unique time steps
taus = sorted(df["tau"].unique())


# Time evolution plot

plt.figure(figsize=(12, 5))

for tau in taus:
    sub = df[df["tau"] == tau]
    plt.plot(sub["t"], sub["x"], label=f"x(t), tau={tau}")
    plt.plot(sub["t"], sub["v"], "--", label=f"v(t), tau={tau}")

plt.xlabel("t")
plt.ylabel("x(t) and v(t)")
plt.title("Time Evolution (Improved Euler)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# Phase plot

plt.figure(figsize=(6, 6))

for tau in taus:
    sub = df[df["tau"] == tau]
    plt.plot(sub["x"], sub["v"], label=f"tau={tau}")

plt.xlabel("x")
plt.ylabel("v")
plt.title("Phase Plot (Improved Euler)")
plt.legend()
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.show()
