import numpy as np
import matplotlib.pyplot as plt
import os

steps = [20, 50, 100, 200]

base_path = "../compare_results"

for N in steps:
    exp_file = f"{base_path}/explicit_{N}.txt"
    imp_file = f"{base_path}/implicit_{N}.txt"
    impv_file = f"{base_path}/improved_{N}.txt"

    data_exp = np.loadtxt(exp_file, usecols=(0, 1, 2))
    data_imp = np.loadtxt(imp_file, usecols=(0, 1, 2))
    data_impv = np.loadtxt(impv_file, usecols=(0, 1, 2))

    plt.figure(figsize=(10,5))
    plt.plot(data_exp[:,0], data_exp[:,1], label='Explicit Euler - position')
    plt.plot(data_imp[:,0], data_imp[:,1], label='Implicit Euler - position')
    plt.plot(data_impv[:,0], data_impv[:,1], label='Improved Euler - position')

    plt.xlabel('time')
    plt.ylabel('position')
    plt.title(f'Mass-Spring Position Over Time (N={N})')
    plt.legend()
    plt.grid()

    plt.savefig(f"{base_path}/time_plot_{N}.png")
    plt.close()

    plt.figure(figsize=(10,5))
    plt.plot(data_exp[:,1], data_exp[:,2], label='Explicit Euler')
    plt.plot(data_imp[:,1], data_imp[:,2], label='Implicit Euler')
    plt.plot(data_impv[:,1], data_impv[:,2], label='Improved Euler')

    plt.xlabel('position')
    plt.ylabel('velocity')
    plt.title(f'Phase Plot (N={N})')
    plt.legend()
    plt.grid()

    plt.savefig(f"{base_path}/phase_plot_{N}.png")
    plt.close()




