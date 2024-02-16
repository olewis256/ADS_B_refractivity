import numpy as np
import matplotlib.pyplot as plt
import smplotlib

data = np.genfromtxt('raypath/PAPERII_paths1.txt')

s = data[:,0]
hf = data[:,1]

sb = data[:,3]
hb = data[:,4]

plt.plot(s, np.array(hf) - np.array(hb[::-1]))
#plt.plot(s, hf)
#plt.plot(s, hb[::-1])
plt.ylabel("Diff (km)")
plt.xlabel("Distance (km)")
plt.title("RK1 diff")
plt.savefig("test_plots/RK1_trace_diff.jpeg")
#print(hp-hpb[::-1])