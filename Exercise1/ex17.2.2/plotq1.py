import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# input data
df = pd.read_csv("explicit_euler.csv")

# unique step sizes
taus = sorted(df["tau"].unique())


# Plot1 : Time evolution plot

plt.figure(figsize=(10, 5))

for tau in taus:
    sub = df[df["tau"] == tau]
    plt.plot(sub["t"], sub["x"], label=f"x(t), tau={tau}")
    plt.plot(sub["t"], sub["v"], label=f"v(t), tau={tau}", linestyle="--")

plt.xlabel("t")
plt.ylabel("x(t) and v(t)")
plt.title("Time Evolution (Explicit Euler)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot 2 : phase plot
plt.figure(figsize=(6, 6))

for tau in taus:
    sub = df[df["tau"] == tau]
    plt.plot(sub["x"], sub["v"], label=f"tau={tau}")

plt.xlabel("x")
plt.ylabel("v")
plt.title("Phase Plot (Explicit Euler)")
plt.legend()
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.show()
