import pandas as pd
import matplotlib.pyplot as plt

df_explicit = pd.read_csv("xian.csv")
df_implicit = pd.read_csv("yin.csv")
df_cn = pd.read_csv("cn.csv")

value_col = "UC" if "UC" in df_explicit.columns else "y"
taus = sorted(df_explicit["tau"].unique())

for tau in taus:
    plt.figure(figsize=(10, 5))

    df_e = df_explicit[df_explicit["tau"] == tau]
    df_i = df_implicit[df_implicit["tau"] == tau]
    df_c = df_cn[df_cn["tau"] == tau]

    plt.plot(df_e["t"], df_e[value_col], label=f"Explicit, τ={tau}", linewidth=1.5)
    plt.plot(df_i["t"], df_i[value_col], label=f"Implicit, τ={tau}", linewidth=1.5)
    plt.plot(df_c["t"], df_c[value_col], label=f"CN, τ={tau}", linewidth=1.5)

    plt.title(f"Comparison for tau = {tau}")
    plt.xlabel("t")
    plt.ylabel(value_col)
    plt.grid(True)
    plt.legend()
    plt.show()
