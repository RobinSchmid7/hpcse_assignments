"""Author: Robin Schmid, 16-927-725, schmirob@ethz.ch"""

import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('transpose_times.txt', comments='#', delimiter=' ', skiprows=1)
block_sizes = [2, 4, 8, 16, 32, 64, 128]

# Selecting unoptimized times
unopt_1 = np.zeros(7)
unopt_2 = np.zeros(7)
unopt_3 = np.zeros(7)
for i in range(0,7):
    unopt_1[i] = data[0][1]
    unopt_2[i] = data[1][1]
    unopt_3[i] = data[2][1]

# Plotting for unoptimized
plt.plot(block_sizes, unopt_1, 'r--', label='N = 1024, unoptimized')
plt.plot(block_sizes, unopt_2, 'g--', label='N = 2048, unoptimized')
plt.plot(block_sizes, unopt_3, 'b--', label='N = 4096, unoptimized')
# Plotting for N = 1024
plt.plot(block_sizes, data[0][2:9], 'r', label='N = 1024')
# Plotting for N = 2048
plt.plot(block_sizes, data[1][2:9], 'g', label='N = 2048')
# Plotting for N = 4096
plt.plot(block_sizes, data[2][2:9], 'b', label='N = 4096')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Block Size')
plt.ylabel('Time [s]')
plt.xticks(block_sizes, block_sizes)
plt.legend()
plt.title('Transposition with different block sizes')
plt.grid(True)
plt.savefig('transpose_plt.pdf')
plt.show()
