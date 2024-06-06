import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import smplotlib
import numpy as np
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
import iris
import iris.plot as iplt
import iris.quickplot as qplt
from scipy.interpolate import UnivariateSpline
import argparse
import os

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


data = iris.load("UKV/extracted_fields_oct25_south.nc")
h = iris.load_cube("UKV/orography_south.nc")

cube0 = data[0][:,0:70,:,:]
cube1 = data[1][:,0:70,:,:]
cube2 = data[2][:,0:70,:,:]

P = cube0.data/100
T = cube1.data
q = cube2.data

e = q*P / ( 0.622 + q - q*0.622)

Nd = 77.6*P/T
Nw = 3.73e5*e/(T**2)

N = Nd + Nw

lat_Clee = 52.39687345348266
lon_Clee = -2.594612181110821

# xc, yc = np.meshgrid([lon_Clee], [lat_Clee])

# lonc, latc = iris.analysis.cartography.rotate_pole(xc, yc, 177.5, 37.5)
# xc, yc = iris.analysis.cartography.unrotate_pole(lonc, latc, 177.5, 37.5)

# print(xc, yc)

glat = cube0.coord('grid_latitude').points
glon = cube0.coord('grid_longitude').points

# x, y = np.meshgrid(glon, glat)

# x, y  = iris.analysis.cartography.unrotate_pole(x, y, 177.5, 37.5)


hlat = h.coord('grid_latitude').points
hlon = h.coord('grid_longitude').points

x, y = np.meshgrid(glon, glat)
xh, yh = np.meshgrid(hlon, hlat)

lons, lats = iris.analysis.cartography.unrotate_pole(x, y, 177.5, 37.5)
lonsh, latsh = iris.analysis.cartography.unrotate_pole(xh, yh, 177.5, 37.5)

point = (lon_Clee, lat_Clee)

dist = np.sqrt((lons - point[0])**2 + (lats - point[1])**2)

min_index = np.unravel_index(np.argmin(dist), dist.shape)

xlon = min_index[1]
xlat = min_index[0]

parser = argparse.ArgumentParser(
                    prog='Generate refractivity data from UKV',
                    description='Tools to generate, plot and save refractivity data generated from UKV',
                    epilog='')

parser.add_argument('--profile', action='store_true',
                    help='Generate a 1D vertical refractivity profile')

parser.add_argument('--profile_grad', action='store_true',
                    help='Generate a 1D vertical refractivity gradient profile')

parser.add_argument('--plot', action='store_true',
                    help='Plot the data')

parser.add_argument('--time',type=int,
                    help='Choose time (UTC)')

args = parser.parse_args()

t = args.time

# ax = plt.axes(projection=ccrs.PlateCarree())

# # define the coordinate system that the grid lons and grid lats are on

# plot = plt.pcolormesh(lons, lats, h.data, cmap='rainbow')
# plt.scatter(lon_Clee, lat_Clee, s=10)
# plt.colorbar(plot)
# ax.coastlines()

# plt.show()


if(args.profile):

    for i in range(62):
        
        cube0.coord('level_height').points[i] += h[xlat][xlon].data * (1.0 - (cube0.coord('level_height').points[i]/40e3) / 0.4338236) ** 2

    if(args.plot):
        plt.plot(N[t,:,xlat,xlon], cube0.coord('level_height').points)
        # plt.plot(Nd[t,:,xlat,xlon], cube0.coord('level_height').points)
        plt.ylim(0, 12e3)
        plt.show()

        

    full_array = np.stack([cube0.coord('level_height').points, N[t,:,xlat,xlon], Nd[t,:,xlat,xlon], Nw[t,:,xlat,xlon]], axis=1)

    path = "../refractivity/OCT25_N_{}.00_(km)_updated_orog.txt".format(t)

    if(os.path.exists(path)):
    
        ftrue = input("Data already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving refractivity file to: {}".format(path))

        np.savetxt(path, full_array, delimiter=" ")

if(args.profile_grad):

    for i in range(62):
        
        cube0.coord('level_height').points[i] += h[xlat][xlon].data * (1.0 - (cube0.coord('level_height').points[i]/40e3) / 0.4338236) ** 2

    if(args.plot):
        plt.plot(q[t,:,xlat,xlon], cube0.coord('level_height').points)
        # plt.plot(Nd[t,:,xlat,xlon], cube0.coord('level_height').points)
        plt.ylim(0, 12e3)
        plt.show()

    interp_func = UnivariateSpline(cube0.coord('level_height').points, N[t,:,xlat,xlon])
    grad = interp_func.derivative()


    full_array = np.stack([cube0.coord('level_height').points, grad(cube0.coord('level_height').points)], axis=1)


    path = "../refractivity/OCT25_N_{}.00_(km)_updated_orog_grad.txt".format(t)

    if(os.path.exists(path)):
    
        ftrue = input("Data already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving refractivity gradient file to: {}".format(path))

        np.savetxt(path, full_array, delimiter=" ")


# ax = fig.add_subplot()

# #Ngrad = interp_func(d_interp, cube0.coord('level_height').points)

# else:

#     direction = 45

#     lonp = [lons[xlat][xlon]]
#     latp = [lats[xlat][xlon]]

#     cubep0 = cube0.copy()

#     Ngrad = np.zeros(70)
    

#     for k in range(62):
                    
#         cubep0.coord('level_height').points.data[k] += h[xlat][xlon].data * (1.0 - (cubep0.coord('level_height').points.data[k]/40e3) / 0.4338236) ** 2
    
#     Ngrad[0] = (N[t,1,xlat,xlon] - N[t,0,xlat,xlon]) / (cubep0.coord('level_height').points.data[1] - cubep0.coord('level_height').points.data[0])
                 
#     for k in range(1, 69):
#         Ngrad[k] = (N[t,k-1,xlat,xlon] - N[t,k+1,xlat,xlon]) / (cubep0.coord('level_height').points.data[k-1] - cubep0.coord('level_height').points.data[k+1])
            
#     Ngrad[69] = (N[t,69,xlat,xlon] - N[t,68,xlat,xlon]) / (cubep0.coord('level_height').points.data[69] - cubep0.coord('level_height').points.data[68])
                           
             

#     hp = [cubep0.coord('level_height').points.data]
#     dp = [0]
#     hs = [h[xlat][xlon].data]
#     Npgrad = [Ngrad]

#     Np, Nwp, Ndp = [N[t,:,xlat,xlon].data], [Nw[t,:,xlat,xlon].data], [Nd[t,:,xlat,xlon].data]


    
#     for i in range(320):
#         for j in range(333):

#             if direction+2 > Bearing(lons[j,i], lats[j,i]) > direction-2:
                
#                 lonp.append(lons[j,i])
#                 latp.append(lats[j,i])

#                 cubep = data[0][:,0:70,:,:]

#                 for k in range(62):
                    
#                     cubep.coord('level_height').points.data[k] += h[j][i].data * (1.0 - (cubep.coord('level_height').points.data[k]/40e3) / 0.4338236) ** 2
                
#                 Ngrad[0] = (N[t,1,j,i] - N[t,0,j,i]) / (cubep.coord('level_height').points.data[1] - cubep.coord('level_height').points.data[0])
                 
#                 for k in range(1, 69):
#                     Ngrad[k] = (N[t,k-1,j,i] - N[t,k+1,j,i]) / (cubep.coord('level_height').points.data[k-1] - cubep.coord('level_height').points.data[k+1])
            
#                 Ngrad[69] = (N[t,69,j,i] - N[t,68,j,i]) / (cubep.coord('level_height').points.data[69] - cubep.coord('level_height').points.data[68])
                 

#                 hp.append(cubep.coord('level_height').points.data)
#                 hs.append(h[j][i].data)
#                 dp.append(Distance(lons[j,i], lats[j,i]))
                
#                 Np.append(N[t,:,j,i].data)
#                 Nwp.append(Nw[t,:,j,i].data)
#                 Ndp.append(Nd[t,:,j,i].data)

#                 Npgrad.append(Ngrad)
                
#                 break


# X = np.tile(np.array(dp), (np.array(hp).shape[1], 1)).T

# print(np.shape(dp))
# print(np.shape(hp))
# print(np.shape(Np))

# fig = plt.figure(figsize=[8, 12])

# ax1 = fig.add_subplot(projection=ccrs.PlateCarree())
# lines = ax1.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False,linewidth=0.1)

# lines.xlabels_top = False
# lines.ylabels_right = False
# ax1.coastlines(zorder=3)
# co = ax1.contourf(lons, lats, h.data,
#                 transform=ccrs.PlateCarree(),
#                 levels=50, cmap='rainbow')
# ax1.plot(lonp, latp, transform=ccrs.PlateCarree())
# ax1.scatter(lon_Clee, lat_Clee, transform=ccrs.PlateCarree(), marker='*', s=80)

# cbar = fig.colorbar(co,ax=ax1, orientation='horizontal', fraction=.08, pad=0.04,aspect=12)
# cbar.set_label('Elevation (m)', labelpad=0.2)

# ax1.set_extent([ -6,  5, 49,  60], crs=ccrs.PlateCarree())
# plt.tight_layout()
# plt.savefig("../plots/surface_topog.jpeg", dpi=700)
# print("test")
# plt.show()

# fig = plt.figure(figsize=[10, 8])

# ax1 = fig.add_subplot(131, projection=ccrs.PlateCarree())

# ax1.coastlines(zorder=3)
# co = ax1.contourf(lons, lats, h.data,
#                 transform=ccrs.PlateCarree(),
#                 levels=50, alpha=0.2, vmin =300)
# ax1.scatter(lonp, latp, transform=ccrs.PlateCarree())
# ax1.scatter(-2.7139, 52.3677, transform=ccrs.PlateCarree())
# cbar = fig.colorbar(co,ax=ax1, orientation='horizontal', fraction=.08, pad=0.04, shrink=0.8, aspect=12)
# cbar.set_label('Elevation (km)', labelpad=0.2)

# ax1.set_aspect('equal')

# ax2 = fig.add_subplot(132)

# ax2.plot(dp, hs)
# ax2.set_xlabel("Distance (km)")
# ax2.set_ylabel("Altitude (m)")

# ax3 = fig.add_subplot(133)
# fig = plt.figure(figsize=[40, 20])

# hp =np.array(hp)/1e3

# d_interp = np.linspace(0, 247, 247)


# print("About to interpolate...")
# #interp_func = interp2d(d_interp,cube0.coord('level_height').points, np.array(Np).T, kind='cubic')

# ax = fig.add_subplot()

# #Ngrad = interp_func(d_interp, cube0.coord('level_height').points)
             
# plot = ax.contourf(X, hp, np.array(Npgrad)*1e3, levels=300, cmap='rainbow')
# ax.scatter(np.full((70,), 50), hp[0,:])
# ax.set_ylabel("Altitude (km)")
# ax.set_xlabel("Distance")  
# cbar = plt.colorbar(plot)
# cbar.set_label('Horizontal refractivity gradient (ppm/km)', labelpad=20)  
# # ax.set_ylim(0,12)

# plt.savefig("../plots/N2D_grad.jpeg", dpi=700)
    
# plt.show()
