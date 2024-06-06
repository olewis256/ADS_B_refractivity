import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import imageio
import math
import smplotlib
import cartopy.crs as ccrs

lat_Clee = 52.39687345348266
lon_Clee = -2.594612181110821

r0 = 6371.229

def Bearing(lon_s, lat_s):
    """ returns bearing direction """
    
    lon_s = lon_s*np.pi/180
    lat_s = lat_s*np.pi/180
    
    X = np.cos(lat_s) * np.sin(lon_s - lon_Clee*np.pi/180)
    Y = (np.cos(lat_Clee*np.pi/180) * np.sin(lat_s)) - (np.sin(lat_Clee*np.pi/180) * np.cos(lat_s) * np.cos(lon_s - lon_Clee*np.pi/180))
    a = np.mod(((math.atan2(X, Y)*180/np.pi)+360), 360) 
  
    return a

def Distance(lon_s, lat_s):
    
    lon_s = lon_s*np.pi/180
    lat_s = lat_s*np.pi/180
    
    A = np.sin((lat_s-lat_Clee*np.pi/180)/2)**2 
    + ( np.cos(lat_Clee*np.pi/180) * np.cos(lat_s) *
       np.sin((lon_s-lon_Clee*np.pi/180)/2)**2 )
    
    d = 2*r0 *np.arcsin(np.sqrt(A))
    
    return d
    
#----------------------------------------------------------------
#                          UKV data
#----------------------------------------------------------------

# max_height_level = 70

# data = Dataset("UKV/extracted_fields_sep21_north-east.nc")
# spec_humid = data.variables['specific_humidity'][:][:,0:max_height_level,:,:]
# air_temp = data.variables['air_temperature'][:][:,0:max_height_level,:,:]
# air_pres = data.variables['air_pressure'][:][:,0:max_height_level,:,:]/100
# time = data.variables['time'][:]
# height = data.variables['level_height'][:]

# e = spec_humid*air_pres / ( 0.622 + 0.378*spec_humid)

# Nd = 77.6*air_pres/air_temp 
# Nw = 3.73e5*e/air_temp**2

# N = Nd + Nw

# t = 16

#----------------------------------------------------------------
#                      Orography data                                  
#----------------------------------------------------------------

datao = Dataset("UKV/orography_sep21_north-east.nc")

lato = datao.variables['grid_latitude'][:]
lono = datao.variables['grid_longitude'][:]
alto = datao.variables['surface_altitude'][:]

lato_copy = lato.copy()
lono_copy = lono.copy()

x = np.zeros((len(lato), len(lono)))
y = np.zeros((len(lato), len(lono)))
z = np.zeros((len(lato), len(lono)))

for i in range(len(lono)):
    for j in range(len(lato)):
        x[j,i] = np.cos(lono[i]*np.pi/180)*np.cos(lato[j]*np.pi/180)
        y[j,i] = np.sin(lono[i]*np.pi/180)*np.cos(lato[j]*np.pi/180)
        z[j,i] = np.sin(lato[j]*np.pi/180)


x_clee = np.cos(lon_Clee*np.pi/180)*np.cos(lat_Clee*np.pi/180)
y_clee = np.sin(lon_Clee*np.pi/180)*np.cos(lat_Clee*np.pi/180)
z_clee = np.sin(lat_Clee*np.pi/180)        
        
rot_lat =  (37.5)*np.pi/180
rot_lon = (180+177.5)*np.pi/180

xn = np.cos(rot_lat)*np.cos(rot_lon)*x + np.sin(rot_lon)*y + np.sin(rot_lat)*np.cos(rot_lon)*z
yn = -np.cos(rot_lat)*np.sin(rot_lon)*x + np.cos(rot_lon)*y - np.sin(rot_lat)*np.sin(rot_lon)*z
zn = -np.sin(rot_lat)*x + np.cos(rot_lat)*z

xn_clee = np.cos(-rot_lat)*np.cos(-rot_lon)*x_clee + np.cos(-rot_lat)*np.sin(-rot_lon)*y_clee + np.sin(-rot_lat)*z_clee
yn_clee = -np.sin(-rot_lon)*x_clee + np.cos(-rot_lon)*y_clee
zn_clee = -np.sin(-rot_lat)*np.cos(-rot_lon)*x_clee - np.sin(-rot_lat)*np.sin(-rot_lon) + np.cos(-rot_lat)*z_clee

lato = np.arcsin(zn)*180/np.pi
lono = np.arctan2(yn, xn)*180/np.pi

lat_clee = np.arcsin(zn_clee)*180/np.pi
lon_clee = np.arctan2(yn_clee, xn_clee)*180/np.pi

lato_arr = lato.reshape(-1, 1)  
lono_arr = lono.reshape(-1, 1)

obs_coords_lon = np.abs(lono_copy - lon_clee).argmin()
obs_coords_lat = np.where(lato_copy == lat_clee)

print(lono_copy, lon_clee, lono)

#-----------------------------------------------------------
#                    Cross-sections
#-----------------------------------------------------------

N_track = []
Nw_track = []
Nd_track = []

N_track.append(N[t,:,obs_coords_lat,obs_coords_lon])
Nw_track.append(Nw[t,:,obs_coords_lat,obs_coords_lon])
Nd_track.append(Nd[t,:,obs_coords_lat,obs_coords_lon])

lon_track = []
lat_track = []
height_ground_track = []

height_track = []
distance_track = []

lon_track.append(lono[obs_coords_lat,obs_coords_lon])
lat_track.append(lato[obs_coords_lat,obs_coords_lon])
height_ground_track.append(alto[obs_coords_lat,obs_coords_lon])

distance_track.append(0)

height_orog_track = np.array(height)
for k in range(62):
    
    height_orog_track[k] = height_orog_track[k] + alto[obs_coords_lat,obs_coords_lon]* (1.0 -
                                    (height_orog_track[k]/40e3) / 0.4338236) ** 2
    
    
height_track.append(height_orog_track[0:max_height_level])

direction = 45

for i in range(200):
    for j in range(279):

        if direction+0.1> Bearing(lono[j,i], lato[j,i]) > direction-0.1:
            
            lon_track.append(lono[j,i])
            lat_track.append(lato[j,i])
            height_ground_track.append(alto[j, i])
            height_orog_track = np.array(height)
            for k in range(62):
                
                height_orog_track[k] = height_orog_track[k] + alto[j,i]* (1.0 -
                                                (height_orog_track[k]/40e3) / 0.4338236) ** 2
                
                
            height_track.append(height_orog_track[0:max_height_level])
            distance_track.append(Distance(lono[j,i], lato[j,i]))
            
            N_track.append(N[t,:,j,i])
            Nw_track.append(Nw[t,:,j,i])
            Nd_track.append(Nd[t,:,j,i])
            
            break
            
            
print(height_track[0])          

fig = plt.figure(figsize=[20, 20])

ax1 = fig.add_subplot(131, projection=ccrs.Mercator())

ax1.coastlines(zorder=3)
co = ax1.contourf(lono, lato, alto,
                transform=ccrs.PlateCarree(),
                cmap='nipy_spectral', levels=50, alpha=0.2)

cbar = fig.colorbar(co,ax=ax1, orientation='horizontal', fraction=.08, pad=0.04, shrink=0.8, aspect=12)
cbar.set_label(r'Elevation (km)', labelpad=0.2)
ax1.scatter(lono[obs_coords_lat, obs_coords_lon], lato[obs_coords_lat, obs_coords_lon], c='w',s=50,transform=ccrs.PlateCarree(), marker='*')
ax1.scatter(lon_track, lat_track, c='black',s=50,transform=ccrs.PlateCarree(), marker='*')
ax1.set_aspect('equal')

ax2 = fig.add_subplot(132)

ax2.plot(distance_track, height_ground_track)
ax2.set_xlabel("Distance (km)")
ax2.set_ylabel("Altitude (m)")

ax3 = fig.add_subplot(133)

D, H = np.meshgrid(distance_track, height[0:max_height_level])             
#print(height + alto[271, 47])                
plot = ax3.contourf(D.T, np.array(height_track)/1e3, N_track, levels=100)
plt.ylabel("Altitude (km)")
plt.xlabel("Distance")  

cbar = plt.colorbar(plot)
cbar.set_label('Refractivity (ppm)', labelpad=20)  
     
#plt.savefig("Refractivity_corrected_18hr.jpeg", dpi=700)

N_track = np.array(N_track)
Nw_track = np.array(Nw_track)
Nd_track = np.array(Nd_track)

D = np.array(D)
H = np.array(H)

N_flatten = N_track.flatten()
Nw_flatten = Nw_track.flatten()
Nd_flatten = Nd_track.flatten()

d_flatten = D.flatten()
h_flatten = H.flatten()

plt.show()
full_array = np.stack([height_track[0], N[t,:,obs_coords_lat,obs_coords_lon], Nd[t,:,obs_coords_lat,obs_coords_lon], Nw[t,:,obs_coords_lat,obs_coords_lon]], axis=1)
#plt.close("all")       
#print(height_track[0]/1e3)      
np.savetxt("../Refractivity_data/Sep21_N_16.00_(km)_updated_orog.txt", full_array, delimiter=" ")
#plt.plot(height_track[0]/1e3, N[t,:,271,47])
#plt.plot(height_track[0]/1e3, Nd[t,:,271,47])
#plt.xlim(0,10)
df = pd.DataFrame({'N': N_flatten, 'Nw': Nd_flatten, 'Nd': Nw_flatten, 'd': d_flatten, 'h': h_flatten/1e3})
#df.to_csv(r'C:/Users/osl202/OneDrive - University of Exeter/PhD/Paper II/.vscode/Refractivity_data/OCT25_UKV_updated.txt', header=None, index=None, sep=' ', mode='a')

