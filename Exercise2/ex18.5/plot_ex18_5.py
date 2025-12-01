import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("ex18_5.csv")

plt.figure(figsize=(10,5))
plt.plot(df["t"], df["theta"], label="theta(t)")
plt.plot(df["t"], df["omega"], label="omega(t)")
plt.xlabel("time (s)")
plt.ylabel("state")
plt.legend()
plt.grid(True)
plt.title("Pendulum State Evolution")
plt.savefig("ex18_5_plot.png", dpi=200)
plt.show()