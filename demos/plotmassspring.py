import numpy as np
import matplotlib.pyplot as plt

data_exp = np.loadtxt('compare_results/explicit_100.txt', usecols=(0, 1, 2))
data_imp = np.loadtxt('compare_results/implicit_100.txt', usecols=(0, 1, 2))
data_impv = np.loadtxt('compare_results/improved_100.txt', usecols=(0, 1, 2))


plt.figure(figsize=(10,5))
plt.plot(data_exp[:,0], data_exp[:,1], label='Explicit Euler - position')
plt.plot(data_imp[:,0], data_imp[:,1], label='Implicit Euler - position')
plt.plot(data_impv[:,0], data_impv[:,1], label='Improved Euler - position')

plt.xlabel('time')
plt.ylabel('position')
plt.title('Mass-Spring System Position Over Time')
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10,5))
plt.plot(data_exp[:,1], data_exp[:,2], label='Explicit Euler')
plt.plot(data_imp[:,1], data_imp[:,2], label='Implicit Euler')
plt.plot(data_impv[:,1], data_impv[:,2], label='Improved Euler')

plt.xlabel('position')
plt.ylabel('velocity')
plt.title('Phase Plot: Mass-Spring System')
plt.legend()
plt.grid()
plt.show()


