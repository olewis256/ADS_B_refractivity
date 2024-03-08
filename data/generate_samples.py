import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import imageio
import smplotlib

import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade
import cartopy.util as cutil
import datetime
import matplotlib.patches as patches
from matplotlib import font_manager as mfonts
import argparse
import os

# -----------
#  Constants
#------------

a, b = 6378.137, 6356.752314245

epsilon = np.sqrt(1 - (b/a)**2)

lat_Clee = 52.398423
lon_Clee = -2.595478
h_Clee = 0.577

#------------------------------------------------

def N(lat):

    N = a / np.sqrt(1 - (epsilon*np.sin(lat))**2)

    return N

parser = argparse.ArgumentParser(
                    prog='Generate ADS-B data from CSV file',
                    description='Tools to generate, plot and save ADS-B data for use in the adjoint inversion model.',
                    epilog='')

parser.add_argument('--paper_real', action='store_true',
                    help='Script to generate the sample datasets used in the paper: t = 0, 900, 1800, 2700 (s)')

parser.add_argument('--geom', action='store_true',
                    help='Plot percentage difference between retrieved and synthetic profile')

parser.add_argument('--profile', action='store_true',
                    help='Plot retrieved profile')

parser.add_argument('--loss', action='store_true',
                    help='Plot loss and RMS statistics')

parser.add_argument('--true', action='store_true',
                    help='Plot true profile using raw ADS-B observations')

parser.add_argument('--paths', action='store_true',
                    help='Plot true and retrieved flightpaths')

args = parser.parse_args()

dataClee = pd.read_csv('Data/cleehill.ringway.day2.csv')
df_orig = pd.DataFrame(dataClee, columns=['ICAO','CALLSIGN','ALTITUDE','LONGITUDE',
                                          'LATITUDE','TIMESTAMP','DISTANCE','BHATDOTRHAT',
                                          'BHATDOTRHATV','ARCANG','RP','RSURF','H','MODEL',
                                          'OBS', 'OBSC', 'RMSC'])
df_orig = df_orig.sort_values('TIMESTAMP')

df_orig['TIMESTAMP'] = df_orig['TIMESTAMP']*0.4

obsAoA_min, obsAoA_max = 0.0, 2.0

df_orig = df_orig[(df_orig['OBSC'] >= 0.*np.pi/180) & (df_orig['OBSC'] <= 2*np.pi/180)]

azim_min, azim_max = -2.5, 2.5

df_orig = df_orig[(-df_orig['BHATDOTRHAT'] >= np.sin(azim_min*np.pi/180)) & (-df_orig['BHATDOTRHAT'] <= np.sin(azim_max*np.pi/180))]

if(args.paper_real):

    print("Subsetting within azim: {} to {}\nand obsAoA: {} to {}".format(azim_min, azim_max, obsAoA_min, obsAoA_max))

    for t in [0, 900, 1800, 2700]:

        df = df_orig.copy()

        df = df[(df['TIMESTAMP'] >= int(t)) & (df['TIMESTAMP'] < int(t+900))]

        path = "../ADS_B_data/sep_NE_paperII_input_central5_t{}_{}.txt".format(int(t), int(t+900))

        if(os.path.exists(path)):
        
            ftrue = input("File {} found, do you want to overwrite? [y/n]".format(path))

            if(ftrue == 'y'):
            
                os.remove(path)

                path = "../ADS_B_data/sep_NE_paperII_input_central5_t{}_{}.txt".format(int(t), int(t+900))
            
        else:

            print("File not found, creating new file: {}".format(path))

        obsAoA = np.arcsin(np.sin(df['OBSC']))*180/np.pi
        azim = np.arcsin(df['BHATDOTRHAT'])*180/np.pi
        dist = df['DISTANCE']
        height = df['H']
        repAoA = np.arcsin(df['BHATDOTRHATV'])*180/np.pi
        icao = df['ICAO']
        time = df['TIMESTAMP']
        error = df['RMSC']
        radius = df['RP']
        lons = df['LONGITUDE']
        lats = df['LATITUDE']
        arcang = df['ARCANG']
        rsurf = df['RSURF']

        lon_i, lat_i = lons*np.pi/180, lats*np.pi/180

        Z_offset = np.sin(lat_Clee*np.pi/180)*N(lat_Clee*np.pi/180)*epsilon**2

        # -------------------------------
        #  Observer geodetic coordinates
        #--------------------------------

        Xog = (N(lat_Clee*np.pi/180) + h_Clee)*np.cos(lat_Clee*np.pi/180)*np.cos(lon_Clee*np.pi/180)
        Yog = (N(lat_Clee*np.pi/180) + h_Clee)*np.cos(lat_Clee*np.pi/180)*np.sin(lon_Clee*np.pi/180)
        Zog = (N(lat_Clee*np.pi/180) + h_Clee)*np.sin(lat_Clee*np.pi/180)

        # ------------------------
        #  Geocentric coordinates
        #-------------------------

        X = (N(lat_i) + height)*np.cos(lat_i)*np.cos(lon_i)
        Y = (N(lat_i) + height)*np.cos(lat_i)*np.sin(lon_i)
        Z = (N(lat_i)*(1-epsilon**2) + height)*np.sin(lat_i) 

        # ------------------------------
        #  Observed-centric coordinates
        #-------------------------------

        Xp = X
        Yp = Y
        Zp = (N(lat_i)*(1-epsilon**2) + height)*np.sin(lat_i) + Z_offset

        Rp = np.sqrt(Xp**2 + Yp**2 + Zp**2)
        Rog = np.sqrt(Xog**2 + Yog**2 + Zog**2)

        air_coords = np.array([Xp/Rp, Yp/Rp, Zp/Rp])
        obs_coords = np.array([Xog/Rog, Yog/Rog, Zog/Rog])

        arc_angle_o = np.arccos(np.dot(obs_coords, air_coords))

        dX = Xp - Xog
        dY = Yp - Yog
        dZ = Zp - Zog

        dR = np.sqrt(dX**2 + dY**2 + dZ**2)
        Rog = np.sqrt(Xog**2 + Yog**2 + Zog**2)

        airdir_coords = np.array([dX/dR, dY/dR, dZ/dR])

        obs_coords = np.array([Xog/Rog, Yog/Rog, Zog/Rog])

        zenith = np.arccos(np.dot(obs_coords, airdir_coords))

        cozenith = np.pi/2 - zenith

        df_sample = pd.DataFrame({'obsAoA': obsAoA, 'h': height, 'd': dist, 'repAoA': cozenith*180/np.pi, 'azim': azim, 'timestamp': time, 'arcang': arc_angle_o, 'lat': lats, 'lon': lons, 'icao': icao})

        print("Dataset length: ", len(obsAoA))
        print("N at observer: ", N(lat_Clee*np.pi/180), "km")

        plt.scatter(-azim, obsAoA, color='green',s=3)
        plt.scatter(-azim, repAoA, color='black',s=3)
        plt.scatter(-azim, cozenith*180/np.pi, color='blue',s=3)
        plt.show()

        df_sample.to_csv(path, header=None, index=None, sep=' ', mode='a')

        if(args.geom):

            max_lat = max(lats)/90 * np.pi/2
            min_lat = lat_Clee/90 * np.pi/2
            ang = np.linspace(min_lat, max_lat, 100)

            x = N(ang)*np.cos(ang)
            z = (N(ang)*(1-epsilon**2))*np.sin(ang) + Z_offset
            xo = N(ang)*np.cos(ang)
            zo = (N(ang)*(1-epsilon**2))*np.sin(ang)

            xs = N(lat_Clee*np.pi/180)*np.cos(ang)
            zs = N(lat_Clee*np.pi/180)*np.sin(ang) 

            xsc = N(lat_Clee*np.pi/180)*np.cos(ang)
            zsc = N(lat_Clee*np.pi/180)*(1-epsilon**2)*np.sin(ang)

            plt.plot(ang*180/np.pi, (np.sqrt(xs**2+zs**2)-np.sqrt(x**2+z**2))*1e3, color='green', label='Relative to obs ROC')
            # plt.plot(ang*180/np.pi, (np.sqrt(xsc**2+zsc**2)-np.sqrt(xo**2+zo**2))*1e3, color='black', label='Relative to obs geo radius')
            # plt.plot(ang*180/np.pi, np.sqrt(x**2+z**2), color='green', label='ellipsoid')
            plt.legend()
            plt.ylabel("Deviation from ellipsoid (metres)")
            plt.xlabel("Latitude (deg.)")

            plt.show()



# fig = plt.figure(figsize=[12, 12])

# date = datetime.datetime(2023, 11, 27, 12)

# ax1 = fig.add_subplot(1,2,1, projection=ccrs.PlateCarree())
# lines = ax1.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False,linewidth=0.1)
# #lines.xlines = False
# #lines.ylines = False
# lines.xlabels_top = False
# lines.ylabels_right = False
# #ax1.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)
# #dt = datetime.strptime(aurora['Forecast Time'], '%Y-%m-%dT%H:%M:%SZ')
# ax1.set_extent([-6.1, 3, 49, 60], crs=ccrs.PlateCarree())
# ax1.coastlines(zorder=3)
# plot = ax1.scatter(lons, lats,
#                 transform=ccrs.PlateCarree(),
#                 s=7,edgecolors='none', c=dist,cmap='rainbow')

# ax1.scatter(-2.596505, 52.398238, s=150, marker='*', color='black',transform=ccrs.PlateCarree())
# cbar = fig.colorbar(plot, ax=ax1, orientation='horizontal', fraction=.05, pad=0.03, shrink=1, aspect=11)
# cbar.set_label('Distance (km)',labelpad=10)
# plt.tight_layout()

# plt.savefig("../Plots/PAPERII_map_all_2.jpeg", dpi=700)
# plt.show()
# defining functions for scalebar


