import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('legendre_data.txt')
x = data[:, 0]

for k in range(6):
    col_val = 1 + 2 * k  # P_k is in column 1 + 2*k
    col_deriv = 2 + 2 * k  # P'_k is in columns 2 + 2*k
    
    y_val = data[:, col_val]
    y_deriv = data[:, col_deriv]

    plt.figure(figsize=(8, 5))
    
    plt.plot(x, y_val, label=f'$P_{k}(x)$ (Polynomial)', color='blue', linewidth=2)
    plt.plot(x, y_deriv, label=f"$P'_{k}(x)$ (Derivative)", color='red', linestyle='--')
    
    plt.title(f'Legendre Polynomial Degree {k}')
    plt.xlabel('x')
    plt.ylabel('Value')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.axhline(0, color='black', linewidth=0.8)

    filename = f'legendre_deg{k}.png'
    plt.savefig(filename)
    
    plt.close()