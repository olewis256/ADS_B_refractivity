from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from metpy.plots import SkewT
from metpy.units import pandas_dataframe_to_unit_arrays, units
from metpy.calc import surface_based_cape_cin, cape_cin, dewpoint_from_relative_humidity, parcel_profile
import numpy as np
from siphon.simplewebservice.wyoming import WyomingUpperAir
from meteostat import Stations, Point, Daily
import matplotlib.dates as mdates
#import time
#from urllib.error import HTTPError
import pandas as pd
import smplotlib
from scipy.interpolate import CubicSpline

r0 = 6365.362 

dt = datetime(2023, 5,9,00)
station = '3354'

df = WyomingUpperAir.request_data(dt, station)
data = pandas_dataframe_to_unit_arrays(df)

p = data['pressure'].magnitude
T = data['temperature'].magnitude
Td = data['dewpoint'].magnitude
h = data['height'].magnitude

RH = np.exp( 17.625 * Td / (243.04 + Td) ) / np.exp(17.625 * T / (243.04 + T))

buck1 = 23.036 - T/333.7
buck2 = T/(279.82+T)
wvpsat = 6.1115*np.exp(buck1*buck2)

N = 77.6*p/(T + 273.15) + 3.73e5 * wvpsat*RH/(T + 273.15)**2

Nwet = 3.73e5 * wvpsat*RH/(T + 273.15)**2

Ndry = 77.6 * p / (T + 273.15)

for i in range(len(h)-1):
    if h[i+1] < h[i]:
        
        h[i+1] = h[i+1] + 10000
for i in range(len(h)-1):
    if h[i+1] < h[i]:
        
        h[i+1] = h[i+1] + 10000
for i in range(len(h)-1):
    if h[i+1] < h[i]:
        
        h[i+1] = h[i+1] + 10000

fig, ax = plt.subplots()

ax.plot(RH*100, h/1e3,color='green')
# ax.plot(N, h/1e3, label='total (wet + dry) refractivity', color='black')
#ax.plot(RH, h/1e3, label="dewpoint")
plt.legend()
ax.set_ylim(-0.5, 14)
# ax.set_xlim(50, 390)
ax.set_ylabel("Altitude (km)")
ax.set_xlabel("Relative humidity (%)")
#ax.set_title("Camborne (03808): 26/8 12z")
# plt.savefig("refrac_9Sep.jpeg", dpi=400)
full_array = np.stack([h, N, Ndry, Nwet], axis=1)
mask = ~np.isnan(N) & ~np.isnan(h)

Nt = N[mask]
h_copy = h[mask]
cs = CubicSpline(h_copy, Nt)
print(cs(0.577*1e3))
buck1 = 23.036 - 9/333.7
buck2 = 9/(279.82+9)
wvpsat = 6.1115*np.exp(buck1*buck2)
print(N)
#N = 77.6*980/(9 + 273.15) + 3.73e5 * 0.9*wvpsat/(9 + 273.15)**2
print(T)
print(Td)
print(h)
# np.savetxt("../refractivity/9May_00z_2024_Watnall_profile_RH.txt", full_array, delimiter=" ")
plt.show()
# =============================================================================
# data = np.genfromtxt("Radiosonde_Watnall_28_Apr_0930.txt",skip_header=0)
# 
# full_array = np.stack([data[:,0], data[:,5], data[:,6], data[:,1], data[:,4]], axis=1)
# 
# np.savetxt("28_95_April_full.txt", full_array, delimiter=" ")
# =============================================================================
#ax.invert_yaxis()


# fig = plt.figure()

# # Initiate the skew-T plot type from MetPy class loaded earlier
# skew = SkewT(fig, rotation=45)

# # Plot the data using normal plotting functions, in this case using
# # log scaling in Y, as dictated by the typical meteorological plot
# skew.plot(p, T, 'r')
# skew.plot(p, Td, 'g')
# #skew.plot_barbs(p[::3], u[::3], v[::3], y_clip_radius=0.03)

# # Set some appropriate axes limits for x and y
# skew.ax.set_xlim(-40, 40)
# #skew.ax.set_ylim(1020, 250)
# plt.ticklabel_format(style='plain')
# # Add the relevant special lines to plot throughout the figure
# #skew.plot_dry_adiabats(t0=np.arange(233, 533, 10) * units.K,
# #                       alpha=0.25, color='orangered')
# #skew.plot_moist_adiabats(t0=np.arange(233, 400, 5) * units.K,
# #                         alpha=0.25, color='tab:green')
# #skew.plot_mixing_lines(p=np.arange(1000, 99, -20) * units.hPa,
# #                       linestyle='dotted', color='tab:blue')
# skew.ax.set_xlabel('Temperature (°C)')
# skew.ax.set_ylabel('Pressure (hPa)')
# # Add some descriptive titles
# #plt.title('{} Sounding'.format(station), loc='left')
# #plt.title('Valid Time: {}'.format(dt), loc='right')
# plt.legend()
# plt.savefig("refrac_9Sep.jpeg", dpi=700)


