import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load CSV files
exp_df = pd.read_csv("xian.csv")
imp_df = pd.read_csv("yin.csv")
cn_df  = pd.read_csv("cn.csv")

Ns = sorted(exp_df["N"].unique())

for N in Ns:
    e = exp_df[exp_df["N"] == N]
    i = imp_df[imp_df["N"] == N]
    c = cn_df[cn_df["N"] == N]

 
    plt.figure(figsize=(10,5))
    plt.plot(e["t"], e["UC"], "b", label=f"Explicit Euler (blow up), N={N}")
    plt.title(f"Explicit Euler Only (N={N})")
    plt.xlabel("t")
    plt.ylabel("UC")
    plt.yscale("symlog")  
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

 
    plt.figure(figsize=(10,5))
    plt.plot(i["t"], i["UC"], label=f"Implicit, N={N}")
    plt.plot(c["t"], c["UC"], label=f"CN, N={N}")
    plt.title(f"Implicit vs Crank–Nicolson (N={N})")
    plt.xlabel("t")
    plt.ylabel("UC")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

 
    plt.figure(figsize=(10,5))
    plt.plot(e["t"], np.abs(e["UC"]), label=f"Explicit (|UC|), N={N}")
    plt.plot(i["t"], np.abs(i["UC"]), label=f"Implicit (|UC|), N={N}")
    plt.plot(c["t"], np.abs(c["UC"]), label=f"CN (|UC|), N={N}")

    plt.yscale("log")      
    plt.xlabel("t")
    plt.ylabel("|UC| (log scale)")
    plt.title(f"Log Comparison of All Methods (N={N})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
