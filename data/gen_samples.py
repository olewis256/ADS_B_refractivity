import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import smplotlib

import cartopy.crs as ccrs
import os
import yaml
import argparse

# -----------
#  Constants
#------------

a, b = 6378.137, 6356.752314245

epsilon = np.sqrt(1 - (b/a)**2)

lat_Clee = 52.398423
lon_Clee = -2.595478
h_Clee = 0.575

def N(lat):

    N = a / np.sqrt(1 - (epsilon*np.sin(lat))**2)

    return N

def coord_transform(lat_i, lon_i, height) -> tuple:

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
    # not used?

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

#------------------------------------------------

class Dataset():

    def __init__(self, params):
        params = params['Dataset']
        self.data_path: str = params['data_path']

        self.name: str = params['name']
        self.save_dir: str = params['save_dir']
        self.plot_dir: str = params['plot_dir']

        self.seed: int = params['seed']
        self.size: int = params['size']
        
        self.AoAmin = params['AoAmin']
        self.AoAmax = params['AoAmax']
        self.azimin = params['azimin']
        self.azimax = params['azimax']

        self.columns = params['columns']
        
        print("Reading in: ", self.data_path)
        data = pd.read_csv(self.data_path)
        df = pd.DataFrame(data, columns=self.columns)

        df.sort_values('TIMESTAMP')
        df['TIMESTAMP'] *= 0.4

        df = df[(df['OBSC'] >= self.AoAmin*np.pi/180) 
                & (df['OBSC'] <= self.AoAmax*np.pi/180)]
        
        df = df[(-df['BHATDOTRHAT'] >= np.sin(self.azimin*np.pi/180))
                 & (-df['BHATDOTRHAT'] <= np.sin(self.azimax*np.pi/180))]
        
        self.df_froze = df

        self.size = params['size']
     
    def load(self, time: int, dt: int) -> None:

        self.t = time
        self.dt = dt

        self.df = self.df_froze
           
        self.df = self.df[(self.df['TIMESTAMP'] >= self.t) 
                & (self.df['TIMESTAMP'] < self.t+self.dt)]
        
        obsAoA = np.arcsin(np.sin(self.df['OBSC']))*180/np.pi
        azim = np.arcsin(self.df['BHATDOTRHAT'])*180/np.pi
        dist = self.df['DISTANCE']
        height = self.df['H']
        icao = self.df['ICAO']
        time = self.df['TIMESTAMP']
        lons = self.df['LONGITUDE']
        lats = self.df['LATITUDE']
       
        lon_i, lat_i = lons*np.pi/180, lats*np.pi/180

        cozenith, arc_angle_o, dist = coord_transform(lat_i, lon_i, height)

        self.df_sample = pd.DataFrame({'obsAoA': obsAoA, 'h': height, 'd': dist,
                                       'repAoA': cozenith*180/np.pi, 'azim': azim,
                                       'timestamp': time, 'arcang': arc_angle_o,
                                       'lat': lats, 'lon': lons, 'icao': icao})
        
        self.df_sample = self.df_sample.sample(n=self.size, random_state=self.seed)

        self.len_data = len(self.df_sample)

        print("Dataset length for sample time {} to {}: {}".format(self.t, self.t+self.dt, self.len_data))

    def save(self):

        try:
            self.df_sample
        except NameError:
            print("Cannot find dataframe.  Ensure data is loaded first.")

        self.name_file = self.name + "_azim_{}_{}_time_{}_{}_len_{}_sample{}.txt".format(self.azimin,
                                                         self.azimax,
                                                         self.t,
                                                         self.t+self.dt,
                                                         self.len_data,
                                                         self.seed
                                                         )

        path = os.path.join(self.save_dir, self.name_file)

        if(os.path.exists(path)):

            ftrue = input("File {} found, do you want to overwrite? [y/n]".format(path))

            if(ftrue=='y'):
                os.remove(path)

        else:

            print("File not found, creating new file: {}".format(path))  

        self.df_sample.to_csv(path, header=True, index=None, sep=' ', mode='a')

    def plot_map(self, save: bool=False):

        fig = plt.figure(figsize=[12, 12])

        ax1 = fig.add_subplot(1,2,1, projection=ccrs.PlateCarree())
        lines = ax1.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False,linewidth=0.1)

        lines.xlabels_top = False
        lines.ylabels_right = False

        # ax1.set_extent([min(lons)-1, max(lons)+1, min(lats)-1, max(lats)+1], crs=ccrs.PlateCarree())
        ax1.set_extent([-6, 3, 60, 49], crs=ccrs.PlateCarree())
        ax1.coastlines(zorder=3)

        plot = ax1.scatter(self.df_froze['LONGITUDE'], self.df_froze['LATITUDE'],
                        transform=ccrs.PlateCarree(),
                        s=7,edgecolors='none', c=self.df_froze['DISTANCE'],cmap='viridis')
        
        ax1.scatter(-1.2502, 53.0059, s=150, color='black', transform=ccrs.PlateCarree())
        ax1.scatter(lon_Clee, lat_Clee, s=150, marker='*', color='black',transform=ccrs.PlateCarree())

        cbar = fig.colorbar(plot, ax=ax1, orientation='horizontal', fraction=.05, pad=0.03, shrink=1, aspect=11)
        cbar.set_label('Distance (km)',labelpad=10)

        plt.tight_layout()
        if save:
            save_plot = self.name + "_map.jpeg"
            plot_path = os.path.join(self.plot_dir, save_plot)
            print("Saving map plot to: ", plot_path)
            plt.savefig(plot_path, dpi=700)
        plt.show()
    

if __name__ == '__main__':

    with open('gendata.yaml', 'r') as file:
        config = yaml.safe_load(file)
    parser = argparse.ArgumentParser(
                    prog='Generate ADS-B data from CSV file',
                    description='Tools to generate, plot and save ADS-B data for use in the adjoint inversion model. Arc angles and reported AoAs are calculated relative to the observer\'s geodetic radius vector.',
                    epilog='')

    data = Dataset(config)
    
    print("Using config file with namespace: ", config['Dataset']['name'])

    use_config = input("Use azimuth range {} to {}, and AoA range {} to {} from config file, continue? [y]/n "
                       .format(config['Dataset']['azimin'], config['Dataset']['azimax'], config['Dataset']['AoAmin'], 
                        config['Dataset']['AoAmax']))

    if(use_config == 'y'):

        data.load(0, 900)
        data.save()

        data.load(900, 900)
        data.save()

        data.load(1800, 900)
        data.save()

        data.load(2700, 900)
        data.save()
        

    
 
        
        
    
    



    
