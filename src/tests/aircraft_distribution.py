import numpy as np
import matplotlib.pyplot as plt
import smplotlib

data = np.genfromtxt('../../aircraft_dist/Height_distribution.txt')

d = data[:,0]
h = data[:,1]

plt.scatter(d, h)
plt.show()