import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("freefall.csv")

t = data[:,0]
error = data[:,1]

plt.figure()
plt.plot(t, error)
plt.xlabel("Time")
plt.ylabel("Error")
plt.title("Error vs Time (Free Fall)")
plt.grid()
plt.savefig("error_vs_time.png")
plt.show()