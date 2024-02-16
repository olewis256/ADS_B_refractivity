import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import smplotlib

#------------------
# Distribution data
#-------------------

dist_data0_01 = pd.read_csv("../Noise/PAPERII_noise_NE_RK3_NEW_0.01.txt", sep=' ', header=None)
#dist_data0_05 = pd.read_csv(r'.vscode/Distribution_data/PAPERII_distrib_0.05_NE_RK3.txt', sep=' ', header=None)
# =============================================================================
fig, (ax1, ax2) = plt.subplots(ncols=2,sharey=True,figsize=(10, 6))

ax1.hist(dist_data0_01, bins=100, label="no noise", linestyle='--')

fig.tight_layout()

ax1.set_xlim(-0.03,0.03)
ax2.set_xlim(-0.25,0.25)

ax1.set_xlabel("Observed AoA error (deg.) \n (a)")
ax2.set_xlabel("Observed AoA error (deg.) \n (b)")
ax1.set_ylabel("Number of broadcasts")


ax2.hist(dist_data0_01, bins=100, label="no noise", linestyle='--')

ax2.sharey(ax1)

plt.tight_layout()
plt.show()

#plt.savefig("Plots/Error_distribution_RK3.jpeg", dpi=700)