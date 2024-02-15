import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import smplotlib
import numpy as np
import math

import iris
import iris.plot as iplt
import iris.quickplot as qplt

r0 = 6371

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


data = iris.load("UKV/extracted_fields_sep21_north-east.nc")
h = iris.load_cube("UKV/orography_sep21_north-east.nc")

cube0 = data[0][:,0:70,:,:]
cube1 = data[1][:,0:70,:,:]
cube2 = data[2][:,0:70,:,:]

P = cube0.data/100
T = cube1.data
q = cube2.data

e = q*P / ( 0.622 + 0.378*q)

Nd = 77.6*P/T
Nw = 3.73e5*e/T**2

N = Nd + Nw

t = 16

lat_Clee = 52.39687345348266
lon_Clee = -2.594612181110821

glat = cube0.coord('grid_latitude').points
glon = cube0.coord('grid_longitude').points

x, y = np.meshgrid(glon, glat)

lons, lats = iris.analysis.cartography.unrotate_pole(x, y, 177.5, 37.5)

point = (lon_Clee, lat_Clee)

dist = np.sqrt((lons - point[0])**2 + (lats - point[1])**2)

min_index = np.unravel_index(np.argmin(dist, axis=None), dist.shape)

xlon = min_index[0]
xlat = min_index[1]

print(cube0.coord('level_height').points)

choice = 0

if(choice == 1):

    for i in range(62):
        
        cube0.coord('level_height').points[i] += h[xlat][xlon].data * (1.0 - (cube0.coord('level_height').points[i]/40e3) / 0.4338236) ** 2

    plt.plot(N[t,:,xlat,xlon], cube0.coord('level_height').points)
    plt.show()

    full_array = np.stack([cube0.coord('level_height').points, N[t,:,xlat,xlon], Nd[t,:,xlat,xlon], Nw[t,:,xlat,xlon]], axis=1)

    np.savetxt("../Refractivity_data/Sep21_N_{}.00_(km)_updated_orog.txt".format(t), full_array, delimiter=" ")

else:

    direction = 45

    lonp = [lons[xlon][xlat]]
    latp = [lats[xlon][xlat]]

    cubep0 = cube0.copy()

    for k in range(62):
                    
        cubep0.coord('level_height').points[k] += h[xlat][xlon].data * (1.0 - (cubep0.coord('level_height').points[k]/40e3) / 0.4338236) ** 2
                

    hp = [cubep0.coord('level_height').points]
    dp = [0]
    hs = [h[xlat][xlon].data]

    Np, Nwp, Ndp = [N[t,:,xlat,xlon].data], [Nw[t,:,xlat,xlon].data], [Nd[t,:,xlat,xlon].data]
    
    for i in range(320):
        for j in range(333):

            if direction+0.2> Bearing(lons[j,i], lats[j,i]) > direction-0.2:
                
                lonp.append(lons[j,i])
                latp.append(lats[j,i])

                cubep = data[0][:,0:70,:,:]

                for k in range(62):
                    
                    cubep.coord('level_height').points[k] += h[j][i].data * (1.0 - (cubep.coord('level_height').points[k]/40e3) / 0.4338236) ** 2
                    
                    
                hp.append(cubep.coord('level_height').points.data)
                hs.append(h[j][i].data)
                dp.append(Distance(lons[j,i], lats[j,i]))
                
                Np.append(N[t,:,j,i].data)
                Nwp.append(Nw[t,:,j,i].data)
                Ndp.append(Nd[t,:,j,i].data)
                
                break


X = np.tile(np.array(dp), (np.array(hp).shape[1], 1)).T

print(np.shape(dp))
print(np.shape(hp))


fig = plt.figure(figsize=[20, 20])

ax1 = fig.add_subplot(131, projection=ccrs.Mercator())

ax1.coastlines(zorder=3)
co = ax1.contourf(lons, lats, h.data,
                transform=ccrs.PlateCarree(),
                cmap='rainbow', levels=50, alpha=0.2)
ax1.scatter(lonp, latp, transform=ccrs.PlateCarree())
cbar = fig.colorbar(co,ax=ax1, orientation='horizontal', fraction=.08, pad=0.04, shrink=0.8, aspect=12)
cbar.set_label('Elevation (km)', labelpad=0.2)

ax1.set_aspect('equal')

ax2 = fig.add_subplot(132)

ax2.plot(dp, hs)
ax2.set_xlabel("Distance (km)")
ax2.set_ylabel("Altitude (m)")

ax3 = fig.add_subplot(133)

             
plot = ax3.contourf(X, hp, Np, levels=100)
ax3.set_ylabel("Altitude (km)")
ax3.set_xlabel("Distance")  
ax3.set_ylim(0, 5e3)
cbar = plt.colorbar(plot)
cbar.set_label('Refractivity (ppm)', labelpad=20)  
    
plt.show()
