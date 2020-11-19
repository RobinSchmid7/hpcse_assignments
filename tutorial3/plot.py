import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

data = np.loadtxt("example.txt")

N = data.shape[1]-1
u = data[:,1:N+1]
t = data[:,0]
x = np.linspace(0.0 + 0.5/N, 1.0 - 0.5/N, N)

for i in range(len(t)):
	print(i)
	plt.ylim((-1.1,1.1))
	plt.xlabel("x")
	plt.ylabel("u")
	plt.plot(x,u[i,:])
	plt.title(label = "t="+str(t[i]))
	plt.savefig(str(i).zfill(4)+".png")
	plt.clf()
