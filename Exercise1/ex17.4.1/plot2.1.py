import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#input data
df = pd.read_csv("implicit_euler.csv")

taus = sorted(df["tau"].unique())


# Time evolution plots: y(t) and v(t)

plt.figure(figsize=(12, 5))

for tau in taus:
    sub = df[df["tau"] == tau]
    plt.plot(sub["t"], sub["y"], label=f"y(t), tau={tau}")
    plt.plot(sub["t"], sub["v"], "--", label=f"v(t), tau={tau}")

plt.xlabel("t")
plt.ylabel("y(t), v(t)")
plt.title("Time Evolution (Implicit Euler)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()



# Phase plot: (y, v)

plt.figure(figsize=(6, 6))

for tau in taus:
    sub = df[df["tau"] == tau]
    plt.plot(sub["y"], sub["v"], label=f"tau={tau}")

plt.xlabel("y")
plt.ylabel("v")
plt.title("Phase Plot (Implicit Euler)")
plt.legend()
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.show()
