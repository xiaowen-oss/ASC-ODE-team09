

import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('legendre_data.txt')
x = data[:, 0]
y_val = data[:, 1]
y_deriv = data[:, 2]

plt.figure(figsize=(10, 6))

plt.plot(x, y_val, label=r'$P_5(x)$ (Polynomial)', color='blue', linewidth=2)

plt.plot(x, y_deriv, label=r"$P'_5(x)$ (Derivative)", color='red', linestyle='--')

plt.title('Legendre Polynomial Order 5 and its Derivative')
plt.xlabel('x')
plt.ylabel('Value')
plt.grid(True, alpha=0.3)
plt.legend()
plt.axhline(0, color='black', linewidth=0.8)

plt.savefig('legendre_plot.png')
print("Plot saved as 'legendre_plot.png'")