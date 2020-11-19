"""Author: Robin Schmid, 16-927-725, schmirob@ethz.ch"""

import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('matrix_matrix_times.txt', comments='#', delimiter=' ', skiprows=1)
block_sizes = [2, 4, 8, 16, 32, 64, 128]

# For row major
# Plotting for N = 256
plt.plot(block_sizes, data[0][2:9], 'r*--', label='N = 256, row')
# Plotting for N = 512
plt.plot(block_sizes, data[1][2:9], 'g*--', label='N = 512, row')
# Plotting for N = 1014
plt.plot(block_sizes, data[2][2:9], 'b*--', label='N = 1024, row')
# Plotting for N = 2048
plt.plot(block_sizes, data[3][2:9], 'y*--', label='N = 2048, row')

# For column major
# Plotting for N = 256
plt.plot(block_sizes, data[0][9:16], 'r*-', label='N = 256, column')
# Plotting for N = 512
plt.plot(block_sizes, data[1][9:16], 'g*-', label='N = 512, column')
# Plotting for N = 1014
plt.plot(block_sizes, data[2][9:16], 'b*-', label='N = 1024, column')
# Plotting for N = 2048
plt.plot(block_sizes, data[3][9:16], 'y*-', label='N = 2048, column')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Block Size')
plt.ylabel('Time [s]')
plt.xticks(block_sizes, block_sizes)
plt.legend()
plt.title('Matrix matrix multiplication with different block sizes')
plt.grid(True)
plt.savefig('matrix_matrix_plt.pdf')
plt.show()
