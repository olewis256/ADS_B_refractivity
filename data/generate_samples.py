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
h_Clee = 0.552

#------------------------------------------------

def N(lat):

    N = a / np.sqrt(1 - (epsilon*np.sin(lat))**2)

    return N

def mapData(lons, lats, dist, save=False):

    fig = plt.figure(figsize=[12, 12])

    date = datetime.datetime(2023, 11, 27, 12)

    ax1 = fig.add_subplot(1,2,1, projection=ccrs.PlateCarree())
    lines = ax1.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False,linewidth=0.1)

    lines.xlabels_top = False
    lines.ylabels_right = False

    # ax1.set_extent([min(lons)-1, max(lons)+1, min(lats)-1, max(lats)+1], crs=ccrs.PlateCarree())
    ax1.set_extent([-6, 3, 60, 49], crs=ccrs.PlateCarree())
    ax1.coastlines(zorder=3)
    plot = ax1.scatter(lons, lats,
                    transform=ccrs.PlateCarree(),
                    s=7,edgecolors='none', c=dist,cmap='viridis')
    ax1.scatter(-1.2502, 53.0059, s=150, color='black', transform=ccrs.PlateCarree())
    ax1.scatter(lon_Clee, lat_Clee, s=150, marker='*', color='black',transform=ccrs.PlateCarree())
    cbar = fig.colorbar(plot, ax=ax1, orientation='horizontal', fraction=.05, pad=0.03, shrink=1, aspect=11)
    cbar.set_label('Distance (km)',labelpad=10)
    plt.tight_layout()
    
    plt.savefig("../plots/PAPERII_map_NE.jpeg", dpi=700)
    plt.show()

def coord_transform(lat_i, lon_i, height):

    Z_offset = np.sin(lat_Clee*np.pi/180)*N(lat_Clee*np.pi/180)*epsilon**2

    # -------------------------------
    #  Observer geodetic coordinates
    #--------------------------------

    Xog = (N(lat_Clee*np.pi/180) + h_Clee)*np.cos(lat_Clee*np.pi/180)*np.cos(lon_Clee*np.pi/180)
    Yog = (N(lat_Clee*np.pi/180) + h_Clee)*np.cos(lat_Clee*np.pi/180)*np.sin(lon_Clee*np.pi/180)
    Zog = (N(lat_Clee*np.pi/180) + h_Clee)*np.sin(lat_Clee*np.pi/180)

    # -------------------------------
    #  Aircraft geodetic coordinates
    #--------------------------------

    Xpg = (N(lat_i) + height)*np.cos(lat_i)*np.cos(lon_i)
    Ypg = (N(lat_i) + height)*np.cos(lat_i)*np.sin(lon_i)
    Zpg = (N(lat_i) + height)*np.sin(lat_i)

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

    airdir_coords = np.array([dX/dR, dY/dR, dZ/dR])

    obs_coords = np.array([Xog/Rog, Yog/Rog, Zog/Rog])

    zenith = np.arccos(np.dot(obs_coords, airdir_coords))

    cozenith = np.pi/2 - zenith

    return cozenith, arc_angle_o, dR

parser = argparse.ArgumentParser(
                    prog='Generate ADS-B data from CSV file',
                    description='Tools to generate, plot and save ADS-B data for use in the adjoint inversion model. Arc angles and reported AoAs are calculated relative to the observer\'s geodetic radius vector.',
                    epilog='')

parser.add_argument('--paper_real', action='store_true',
                    help='Script to generate the sample datasets used in the paper: t = 0, 900, 1800, 2700 (s)')

parser.add_argument('--geom', action='store_true',
                    help='Plot the altitude variation between spherical approximation and ellipsoid')

parser.add_argument('--profile', action='store_true',
                    help='Plot retrieved profile')

parser.add_argument('--loss', action='store_true',
                    help='Plot loss and RMS statistics')

parser.add_argument('--true', action='store_true',
                    help='Plot true profile using raw ADS-B observations')

parser.add_argument('--paths', action='store_true',
                    help='Plot true and retrieved flightpaths')

parser.add_argument('--rand_samples', action='store_true',
                    help='Script to generate random samples using September NE data')

parser.add_argument('--south', action='store_true',
                    help='Script to generate south data (IGARSS)')

parser.add_argument('--synthetic', action='store_true',
                    help='Script to generate synthetic subsample data for September NE')

parser.add_argument('--may', action='store_true',
                    help='Script to generate May subsamples.')

parser.add_argument('--time', type=int,
                    help='Start time (s)')

parser.add_argument('--alt', action='store_true',
                    help='Generate Clee Hill altitude statistics')

args = parser.parse_args()

if(args.alt):
    df_orig = pd.read_csv('Data/clee_hill_alt.txt')
    height = df_orig['altHAE']
    print("Mean height ", np.mean(height))
    print("STD ", np.std(height))

if(args.paper_real):

    df_orig = pd.read_csv('Data/cleehill.ringway.day2.csv')
    df_orig = pd.DataFrame(df_orig, columns=['ICAO','CALLSIGN','ALTITUDE','LONGITUDE',
                                            'LATITUDE','TIMESTAMP','DISTANCE','BHATDOTRHAT',
                                            'BHATDOTRHATV','ARCANG','RP','RSURF','H','MODEL',
                                            'OBS', 'OBSC', 'RMSC'])
    df_orig = df_orig.sort_values('TIMESTAMP')

    df_orig['TIMESTAMP'] = df_orig['TIMESTAMP']*0.4
    print(min(df_orig['LATITUDE']), max(df_orig['LATITUDE']), min(df_orig['LONGITUDE']), max(df_orig['LONGITUDE']))
    obsAoA_min, obsAoA_max = 0.0, 2.0

    df_orig = df_orig[(df_orig['OBSC'] >= obsAoA_max*np.pi/180) & (df_orig['OBSC'] <= obsAoA_min*np.pi/180)]

    azim_min, azim_max = -5, 5

    df_orig = df_orig[(-df_orig['BHATDOTRHAT'] >= np.sin(azim_min*np.pi/180)) & (-df_orig['BHATDOTRHAT'] <= np.sin(azim_max*np.pi/180))]

    num = 1

    if(args.rand_samples):

        num = input("Choose number of batches from each data time slot: ")

    num = int(num)

    azim_range = azim_max - azim_min

    print("Subsetting within azim: {} to {}\nand obsAoA: {} to {}".format(azim_min, azim_max, obsAoA_min, obsAoA_max))

    for t in [0, 900, 1800, 2700]:

        df_select = df_orig.copy()

        df_select = df_select[(df_select['TIMESTAMP'] >= int(t)) & (df_select['TIMESTAMP'] < int(t+900))]

        for i in range(num):

            if num == 1:

                path = "../ADS_B_data/sep_NE_paperII_input_central{}_t{}_{}.txt".format(azim_range, int(t), int(t+900))

                df = df_select

            else:

                path = "../ADS_B_data/sep_NE_paperII_input_central{}_t{}_{}_sample_{}.txt".format(azim_range, int(t), int(t+900), i)

                df = df_select.sample(n = int(len(df_select)/num))

                print("Number of ADS-B messages in batch: ", int(len(df_select)/num))

            if(os.path.exists(path)):
            
                ftrue = input("File {} found, do you want to overwrite? [y/n]".format(path))

                if(ftrue == 'y'):
                
                    os.remove(path)
                
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

            cozenith, arc_angle_o, dist = coord_transform(lat_i, lon_i, height)

            df_sample = pd.DataFrame({'obsAoA': obsAoA, 'h': height, 'd': dist, 'repAoA': cozenith*180/np.pi, 'azim': azim, 'timestamp': time, 'arcang': arc_angle_o, 'lat': lats, 'lon': lons, 'icao': icao})

            print("Dataset length: ", len(obsAoA))
            print("N at observer: ", N(lat_Clee*np.pi/180), "km")

            df_sample.to_csv(path, header=None, index=None, sep=' ', mode='a')

            mapData(df_orig['LONGITUDE'], df_orig['LATITUDE'], df_orig['DISTANCE'])

            plt.close()            

            if(args.geom):

                Z_offset = np.sin(lat_Clee*np.pi/180)*N(lat_Clee*np.pi/180)*epsilon**2

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

if(args.south):
    
    df_orig = pd.read_csv('Data/cleesouth.azbase.integrated.csv')
    # df_orig = pd.DataFrame(dataClee, columns=['ICAO','CALLSIGN','ALTITUDE','LONGITUDE',
    #                                         'LATITUDE','TIMESTAMP','DISTANCE','BHATDOTRHAT',
    #                                         'BHATDOTRHATV','ARCANG','RP','RSURF','H','MODEL',
    #                                         'OBS', 'OBSC', 'RMSC'])
    df_orig = df_orig.sort_values('TIMESTAMP')

    df_map = df_orig.copy()

    df_orig['TIMESTAMP'] = df_orig['TIMESTAMP']*0.4

    obsAoA_min, obsAoA_max = 0.0, 2.0

    df_orig = df_orig[(df_orig['OBSC'] >= obsAoA_min*np.pi/180) & (df_orig['OBSC'] <= obsAoA_max*np.pi/180)]

    azim_min, azim_max = -2.5, 2.5

    df_orig = df_orig[(-df_orig['BHATDOTRHAT'] >= np.sin(azim_min*np.pi/180)) & (-df_orig['BHATDOTRHAT'] <= np.sin(azim_max*np.pi/180))]

    print("Subsetting south data within azim: {} to {}\nand obsAoA: {} to {}".format(azim_min, azim_max, obsAoA_min, obsAoA_max))
    
    for t in [0, 900, 1800, 2700]:

        df = df_orig.copy()

        df = df[(df['TIMESTAMP'] >= int(t)) & (df['TIMESTAMP'] < int(t+900))]

        azim_range = azim_max - azim_min

        path = "../ADS_B_data/oct_S_paperII_input_central{}_t{}_{}.txt".format(azim_range, int(t), int(t+900))

        if(os.path.exists(path)):
        
            ftrue = input("File {} found, do you want to overwrite? [y/n]".format(path))

            if(ftrue == 'y'):
            
                os.remove(path)

                path = "../ADS_B_data/oct_S_paperII_input_central5_t{}_{}.txt".format(int(t), int(t+900))
            
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

        cozenith, arc_angle_o, dist = coord_transform(lat_i, lon_i, height)

        df_sample = pd.DataFrame({'obsAoA': obsAoA, 'h': height, 'd': dist, 'repAoA': cozenith*180/np.pi, 'azim': azim, 'timestamp': time, 'arcang': arc_angle_o, 'lat': lats, 'lon': lons, 'icao': icao})

        print("Dataset length: ", len(obsAoA))
        print("N at observer: ", N(lat_Clee*np.pi/180), "km")
        if t ==0:
            plt.scatter(dist, height, c=lats, edgecolors="none") 
            plt.colorbar()
            # plt.scatter(-azim, obsAoA, color='b')
            plt.show()

    mapData(df_map['LONGITUDE'], df_map['LATITUDE'], df_map['DISTANCE'])

    plt.close()

    df_map.sort_values('DISTANCE', ascending=False)

    plot = plt.scatter(np.arcsin(df_map['BHATDOTRHATV'])*180/np.pi, np.arcsin(np.sin(df_map['OBSC']))*180/np.pi, c=df_map['DISTANCE'], edgecolors='none', s=1.5, cmap='viridis')
    cbar = plt.colorbar(plot)
    cbar.set_label("Distance (km)")
    plt.plot([0,2],[0,2], color='black')
    plt.ylim(0.01, 2)
    plt.xlim(0.01, 2)
    plt.xlabel("Reported AoA (deg.)")
    plt.ylabel("Observed AoA (deg.)")
    plt.yscale('log')
    plt.xscale('log')

    plt.savefig("../plots/obs_vs_rep.jpeg", dpi=400)
    plt.show()

        # df_sample.to_csv(path, header=None, index=None, sep=' ', mode='a')


if(args.synthetic):

    dataClee = pd.read_csv('Data/cleehill.ringway.day2.csv')
    df = pd.DataFrame(dataClee, columns=['ICAO','CALLSIGN','ALTITUDE','LONGITUDE',
                                            'LATITUDE','TIMESTAMP','DISTANCE','BHATDOTRHAT',
                                            'BHATDOTRHATV','ARCANG','RP','RSURF','H','MODEL',
                                            'OBS', 'OBSC', 'RMSC'])
    df = df.sort_values('TIMESTAMP')

    df['TIMESTAMP'] = df['TIMESTAMP']*0.4
    print(min(df['LATITUDE']), max(df['LATITUDE']), min(df['LONGITUDE']), max(df['LONGITUDE']))
    obsAoA_min, obsAoA_max = 0.0, 2.0

    df = df[(df['OBSC'] >= obsAoA_min*np.pi/180) & (df['OBSC'] <= obsAoA_max*np.pi/180)]

    azim_min, azim_max = -5, 5

    df = df[(-df['BHATDOTRHAT'] >= np.sin(azim_min*np.pi/180)) & (-df['BHATDOTRHAT'] <= np.sin(azim_max*np.pi/180))]


    print("Subsetting within azim: {} to {}\nand obsAoA: {} to {}".format(azim_min, azim_max, obsAoA_min, obsAoA_max))

    azim_range = azim_max - azim_min

    path = "../ADS_B_data/sep_NE_paperII_input_5.000_central{}.txt".format(azim_range)

        
    if(os.path.exists(path)):
    
        ftrue = input("File {} found, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)
        
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

    cozenith, arc_angle_o, dist = coord_transform(lat_i, lon_i, height)

    df_sample = pd.DataFrame({'obsAoA': obsAoA, 'h': height, 'd': dist, 'repAoA': cozenith*180/np.pi, 'azim': azim, 'timestamp': time, 'arcang': arc_angle_o, 'lat': lats, 'lon': lons, 'icao': icao})

    df_sample = df_sample.sample(n=5000, random_state= 8520)

    print("Dataset length: ", len(df_sample['obsAoA']))
    print("N at observer: ", N(lat_Clee*np.pi/180), "km")

    plt.hist(df_sample['h'], bins=100, color='black', orientation='horizontal')
    plt.ylabel("Altitude (km)")
    plt.xlabel("Number of broadcasts")
    plt.savefig("../plots/PAPERII_aircraft_dist.jpeg", dpi=400)
    df_sample.to_csv(path, header=None, index=None, sep=' ', mode='a')

    mapData(lons, lats, dist)

    if(args.geom):

        max_lat = max(lats)/90 * np.pi/2
        min_lat = lat_Clee/90 * np.pi/2
        ang = np.linspace(min_lat, max_lat, 100)

        x = N(ang)*np.cos(ang)
        z = (N(ang)*(1-epsilon**2))*np.sin(ang) + epsilon**2*N(lat_Clee*np.pi/180)*np.sin(lat_Clee*np.pi/180)
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


if(args.may):

    t = args.time

    df = pd.read_csv('Data/may8data.spacetime.dat', sep="\s+")
    
    df = df.sort_values('TIMESTAMP')

    df['TIMESTAMP'] = df['TIMESTAMP']/5

    print(min(df['LATITUDE']), max(df['LATITUDE']), min(df['LONGITUDE']), max(df['LONGITUDE']))
    print("Maximum time: ", max(df['TIMESTAMP']))
    obsAoA_min, obsAoA_max = 0.0, 3.0

    df = df[(df['OBSC'] >= obsAoA_min*np.pi/180) & (df['OBSC'] <= obsAoA_max*np.pi/180)]


    azim_min, azim_max = -10, 10

    df = df[(-df['BHATDOTRHAT'] >= np.sin(azim_min*np.pi/180)) & (-df['BHATDOTRHAT'] <= np.sin(azim_max*np.pi/180))]
    df = df[(df['TIMESTAMP'] >= int(t)) & (df['TIMESTAMP'] < int(t+1800))]

    print("Subsetting within azim: {} to {}\nand obsAoA: {} to {}".format(azim_min, azim_max, obsAoA_min, obsAoA_max))

    azim_range = azim_max - azim_min

    path = "../ADS_B_data/May8_NE_paperII_input_5.000_azim_{}_to_{}_t{}_to_{}.txt".format(azim_min, azim_max, t, t+1800)
    # path = "../ADS_B_data/May_NE_paperII_input_all.txt"
        
    if(os.path.exists(path)):
    
        ftrue = input("File {} found, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)
        
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

    cozenith, arc_angle_o, dist = coord_transform(lat_i, lon_i, height)

    df_sample = pd.DataFrame({'obsAoA': obsAoA, 'h': height, 'd': dist, 'repAoA': cozenith*180/np.pi, 'azim': azim, 'timestamp': time, 'arcang': arc_angle_o, 'lat': lats, 'lon': lons, 'icao': icao})

    df_sample = df_sample.sample(n=5000, random_state= 8520)

    print("Dataset length: ", len(df_sample['obsAoA']))
    print("N at observer: ", N(lat_Clee*np.pi/180), "km")
   
    df_sample.to_csv(path, header=None, index=None, sep=' ', mode='a')

    plt.scatter(dist, height, c=lats, edgecolors="none") 
    # plt.colorbar()
    # plt.scatter(-azim, obsAoA, c=height)
    plt.show()

    if(args.paths):
        # plt.scatter(azim, obsAoA, s=5, c='blue', edgecolor='none')
        # plt.scatter(azim, repAoA, s=5, c='black',edgecolor='none')

        plt.scatter(lons, lats, s=2)
        plt.show()

    
    
