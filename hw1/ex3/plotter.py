"""
Project: Plotting for HW 1 Exercise 2b)
Maintainer: Robin Schmid
Email: : schmirob@ethz.ch
"""

import numpy as np
import matplotlib.pyplot as plt

M = 100000000
# Load data
data_mode1 = np.loadtxt('data_mode1', comments='#')
data_mode2 = np.loadtxt('data_mode2', comments='#')
data_mode3 = np.loadtxt('data_mode3', comments='#')

# # Generate times mode 1
# t_tot_1 = []
# k = 0
# for i in range(0,len(data_mode1[:,0])):
#     k += data_mode1[i, 3]
#     t_tot_1.append(k)
# # Generate times mode 2
# t_tot_2 = []
# k = 0
# for i in range(0,len(data_mode2[:,0])):
#     k += data_mode2[i, 3]
#     t_tot_2.append(k)
# # Generate times mode 3
# t_tot_3 = []
# k = 0
# for i in range(0,len(data_mode3[:,0])):
#     k += data_mode3[i, 3]
#     t_tot_3.append(k)

# Plot data
plt.plot(data_mode1[:,1], data_mode1[:, 3], 'r*--', label='mode 1')
plt.plot(data_mode2[:,1], data_mode2[:, 3], 'g*--', label='mode 2')
plt.plot(data_mode3[:,1], data_mode3[:, 3], 'b*--', label='mode 3')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('log N')
plt.ylabel('log M/t_tot')
plt.legend()
plt.title('Cache size variants')
plt.grid(True)
plt.savefig('result.pdf')
plt.show()
