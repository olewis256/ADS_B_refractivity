import matplotlib.pyplot as plt
import pandas as pd
import smplotlib

#-------
# Test
#-------

data_1m = pd.read_csv(r'../../Gradient/PAPERII_grad_0.001km.txt', sep=' ',header=None)
data_100m = pd.read_csv(r'../../Gradient/PAPERII_grad_0.1km.txt', sep=' ',header=None)

data_1m.columns = ["FD_grad", "adj_grad", "h", "i"]
data_100m.columns = ["FD_grad", "adj_grad", "h", "i"]

fig, axes = plt.subplots(nrows=2, ncols=2,figsize = (10, 15))
ax1, ax2, ax3, ax4 = axes.flatten()

ax1.plot(abs(data_1m["adj_grad"]), data_1m["h"], linewidth=3,label="adjoint \n(adj)")
ax1.plot(abs(data_1m["FD_grad"]), data_1m["h"], linewidth=3,color='lightgreen',linestyle='--',label="finite difference \n(FD)")
ax1.text(1e3, 9, "(a) 1 m")
ax1.set_ylabel("Altitude (km)")
ax3.set_xlabel(r"$\left|dj/dn\right|$")
ax1.set_ylim(-0.2, 12)

ax1.legend()

ax2.plot(100*abs((data_1m["FD_grad"] - data_1m["adj_grad"]) / (data_1m["FD_grad"] + 1e-13)), data_1m["h"])
ax2.text(2.8, 9, "(b) 1 m")
ax2.set_ylim(-0.2, 12)

ax3.plot(abs(data_100m["adj_grad"]*1), data_100m["h"], linewidth=3)#,label="adjoint \n(adj)")
ax3.plot(abs(data_100m["FD_grad"]), data_100m["h"], linewidth=3,color='lightgreen',linestyle='--')#,label="finite difference \n(FD)")
#ax3.set_xscale("log")
ax3.set_ylabel("Altitude (km)")
ax3.set_xlabel(r"$\left|dj/dn\right|$")
ax3.text(1e3, 9, "(c) 100 m")

ax3.set_ylim(-0.2, 12)

ax4.plot(100*abs((data_100m["FD_grad"] - data_100m["adj_grad"]) / (data_100m["FD_grad"] + 1e-13)), data_100m["h"])
ax4.set_xlabel(r'$ \left|\frac{(dj/dn)_{FD} - (dj/dn)_{adj}}{(dj/dn)_{FD}}\right|$' + ' (%)')
ax4.text(2.8, 9, "(d) 100 m")
ax4.set_ylim(-0.2, 12)
ax1.set_xscale("log")
ax3.set_xscale("log")
ax3.sharex(ax1)
ax2.sharex(ax4)
fig.set_figheight(13)
fig.set_figwidth(9)

plt.tight_layout()

plt.savefig("../../Plots/gradients.jpeg", dpi=700)

