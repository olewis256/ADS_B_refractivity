import numpy as np
import matplotlib.pyplot as plt
import smplotlib

data = np.genfromtxt("../Retrievals/PAPERII_retrieve_NE_RK3_NEW_0.01.txt")

Noptim = data[:,0]
h = data[:,1]
Ntarget = data[:,2]
Ninit = data[:,3]

plt.plot(Noptim, h, color='green', label='retrieved')
plt.plot(Ninit, h, linestyle="--", label='initial')
plt.plot(Ntarget, h, color='black', label='radiosonde')

plt.xlabel("Refractivity")
plt.ylabel("Altitude (km)")

plt.legend()
plt.savefig("../Plots/PAPERII_retrieve_NE_RK3_NEW_0.01_PROFILE.jpeg")
# data2 = np.genfromtxt("PAPERII_paths2.txt")
# data3 = np.genfromtxt("PAPERII_paths1.txt")
# s = data2[:,0][1:3003]
# hp = data2[:,1][1:3003]

# sb = data3[:,0]
# hpb = data3[:,1]

# data_paths = np.genfromtxt("PAPERII_paths1.txt")

# s = data_paths[:,0]
# hf = data_paths[:,1]

# sb = data_paths[:,3]
# hb = data_paths[:,4]

# plt.plot(s, np.array(hf) - np.array(hb[::-1]))
# #plt.plot(s, hf)
# #plt.plot(s, hb[::-1])
# plt.ylabel("Diff (km)")
# plt.xlabel("Distance (km)")
# plt.title("RK1 diff")
# plt.savefig("RK1_trace_diff.jpeg")
# #print(hp-hpb[::-1])


plt.show()