import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('legendre_data.txt')
x = data[:, 0]
max_deg = 5

# This code is to see al polynomials and all derivatives at once
plt.figure(figsize=(10, 6))

colors = ['black', 'blue', 'orange', 'green', 'red', 'purple']

for k in range(max_deg + 1):
    col_val = 1 + 2 * k
    y_val = data[:, col_val]
    plt.plot(x, y_val, label=f'$P_{k}(x)$', color=colors[k], linewidth=2)

plt.title('Legendre Polynomials from degree 0 to 5')
plt.xlabel('x')
plt.ylabel('$P_n(x)$')
plt.grid(True, alpha=0.3)
plt.legend(loc='lower right')
plt.axhline(0, color='black', linewidth=0.8)
plt.ylim(-1.1, 1.1) # Limitamos eje Y para que se vea bien
plt.savefig('legendre_combined_poly.png')

plt.figure(figsize=(10, 6))

for k in range(1, max_deg + 1): # Derivative of P0 is 0 
    col_deriv = 2 + 2 * k
    y_deriv = data[:, col_deriv]
    plt.plot(x, y_deriv, label=f"$P'_{k}(x)$", color=colors[k], linestyle='--')

plt.title('Legendre Polynomials Derivatives (degrees 1 to 5)')
plt.xlabel('x')
plt.ylabel("$P'_n(x)$")
plt.grid(True, alpha=0.3)
plt.legend()
plt.axhline(0, color='black', linewidth=0.8)
plt.savefig('legendre_combined_deriv.png')