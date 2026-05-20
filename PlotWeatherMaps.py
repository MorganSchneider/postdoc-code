# -*- coding: utf-8 -*-
"""
Created on Fri May  1 11:46:48 2026

@author: L. Schielicke
"""

### Import packages ###

from era5utils import *



#%% Choose time, location, region, case (tornado or hail) and file paths

#Didsbury, AB (Tornado)
#Date: 7/1/2023
#Coordinates: (51.6107,-114.198)
#Time (UTC): 17:45

fp = "C:/Users/mschne28/Documents/era5/tor_outbreaks/"

case='tornado'
#case='hail'

# region='Alberta'
#region='SasMan'
region='GreatLakes'
# region='Canada'

yyyyt = 2021
mmt   = 8
ddt   = 11
hht   = 20

# yyyyt = 2021; mmt = 8; ddt = 11; tor hht = 20,21 (18-22)
# yyyyt = 2025; mmt = 6; ddt = 23-24; tor hht = 20,22,1 (19-23, 0-3)
# yyyyt = 2022; mmt = 5; ddt = 30-31; tor hht = 0,1,2 (22-23, 0-3)
# yyyyt = 2022; mmt = 5; ddt = 21; tor hht = 15,17 (13-18)
# yyyyt = 2025; mmt = 7; ddt = 24; bow/db hht = 21-23 (19-23)
# yyyyt = 2026; mmt = 4; ddt = 15; tor/null hht = 5,6 (2-8)

timt="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))


lattstart = 46.4797  #Latitude of tornado/beginning of hail swath (from NTP/NHP database)
lontstart = -83.6403 #Longitude of tornado/beginning of hail swath (from NTP/NHP database)

if (yyyyt==2021) & (mmt==8) & (ddt==11):
    lattstart = [46.4797, 46.5449, 46.5724, 46.4657, 46.3870, 46.6557]
    lontstart = [-83.6403, -83.5907, -83.2889, -83.2155, -82.7493, -82.4617]
elif (yyyyt==2025) & (mmt==6) & ((ddt==23)|(ddt==24)):
    lattstart = [48.3975, 48.4647, 48.2033, 47.9149, 46.7872, 46.7757]
    lontstart = [-75.5880, -75.4900, -73.5234, -73.5056, -70.7827, -70.4072]
elif (yyyyt==2022) & (mmt==5) & ((ddt==30)|(ddt==31)):
    lattstart = [48.6148, 48.6699, 48.9114, 49.3689, 48.6757]
    lontstart = [-93.5304, -93.5219, -93.5778, -93.0754, -92.2308]
elif (yyyyt==2022) & (mmt==5) & (ddt==21):
    lattstart = [43.0179, 42.9217, 44.1058, 44.1755]
    lontstart = [-81.2216, -81.1977, -79.1458, -78.7722]
elif (yyyyt==2025) & (mmt==7) & (ddt==24):
    lattstart = [43.4700, 44.2500, 45.0000]
    lontstart = [-81.1800, -80.5000, -79.2500]
elif (yyyyt==2026) & (mmt==4) & (ddt==15):
    lattstart = [42.2718, 42.2684, 42.2393, 42.2423]
    lontstart = [-83.7524, -83.2104, -83.0576, -82.8484]


latt = 46.50       # latitude of sounding
lont = -83.75     # longitude of sounding #if you do not need this, uncomment the next two lines
latt = np.mean(lattstart)
lont = np.mean(lontstart)

figfolder=fp+'figs/'
data_preslevs=fp+f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_preslevs.nc"
data_singlevs=fp+f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_singlevs.nc"



#% Read data files

# Pressure level data
data = xr.open_dataset(data_preslevs)
data.head(5)

# Single level data
datas = xr.open_dataset(data_singlevs)
datas.head(5)


#%% Make plots

#% Define plot domain - can add other regions here

#Plotdomain:
#lonW=-125
#lonE=-90
#latS=45
#latN=65

# region='Alberta'
#region='SasMan'
region='GreatLakes'
# region='Canada'

if region=='GreatLakes':    
    [lonW,lonE,latS,latN] = [-93.0,-70.0,39.5,52.0]
    n=5
if region=='SasMan':
    [lonW,lonE,latS,latN] = [-114.0,-88.0,46.5,60.5]
    n=5
if region=='Alberta':      
    [lonW,lonE,latS,latN] = [-128.0,-102.0,46.5,60.5]
    n=5
if region=='Canada':
    [lonW,lonE,latS,latN] = [-132.0,-58.0,35.0,75.0]
    n=10



figfolder=None

PlotCAPEMap(data,datas,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,pngfolder=figfolder)

PlotPressureMaps(data,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,pngfolder=figfolder)

PlotMoistureMap2(data,datas,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,pngfolder=figfolder)


data00 = data.sel(valid_time=timt)
z00=data00['z'].values/9.81
u00=data00['u'].values
v00=data00['v'].values
lat=data00['latitude'].values
lon=data00['longitude'].values

# latt=51.5  #Latitude of tornado/beginning of hail swath (from NTP/NHP database)
# lont=-114.0

lati=np.where(lat==latt)
loni=np.where(lon==lont)

# print(np.shape(z00))
# print(lati,loni)

# for i in range(0,np.shape(z00)[0]):
#     print('%d\t%.1f\t%.1f' %(z00[i,lati,loni],u00[i,lati,loni],v00[i,lati,loni]))



#%% Trough axis tilt 2nd attempt - single pressure level

pres = 700

dat1 = data.sel(pressure_level=pres, valid_time=timt)
z1 = dat1["z"].values/9.81
lon1 = dat1["longitude"].values
lat1 = dat1["latitude"].values
u1 = dat1["u"].values
v1 = dat1["v"].values
vort1 = dat1["vo"].values
pv1 = dat1["pv"].values
LON,LAT = np.meshgrid(lon1,lat1)


iyl,ixl = mc.find_peaks(z1, maxima=False)
iyh,ixh = mc.find_peaks(z1)

# dx,dy = mc.lat_lon_grid_deltas(lon1,lat1)

# z1_del2 = mc.laplacian(dat1["z"]/9.81).values

z1_del2 = mc.geospatial_laplacian(dat1["z"]/9.81, crs=ccrs.PlateCarree()).values
z1_del2 = mc.smooth_gaussian(z1_del2, 10)

# z2 = mc.smooth_gaussian(dat1["z"]/9.81, 10)
# z1_del2 = mc.geospatial_laplacian(z2, crs=ccrs.PlateCarree()).values


#%
fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()})
ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linestyle='-')
ax.add_feature(cfeature.STATES, linestyle='-')
ax.add_feature(cfeature.COASTLINE)

if (pres == 850):
    # color_vals = [-40,-32,-28,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,\
    #               8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42]
    xmin_z = 800; xmax_z = 1800; xint_z = 25     # z interval
    color_vals = np.linspace(-1e-8, 1e-8, 41)

elif (pres == 700):
    # color_vals = [-50,-40,-38,-36,-34,-32,-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,\
    #               4,6,8,10,12,14,16,18,20,22,30]
    xmin_z = 2400; xmax_z = 4000; xint_z = 50    # z interval
    color_vals = np.linspace(-1e-8, 1e-8, 41)

elif (pres == 500):
    # color_vals = [-60,-56,-54,-52,-48,-46,-44,-42,-40,-38,-36,-34,-32,-30,-28,-26,-24,-22,\
    #               -20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,10]
    xmin_z = 4800; xmax_z = 6500; xint_z = 50    # z interval
    color_vals = np.linspace(-1e-8, 1e-8, 41)

elif (pres == 300):
    # color_vals = [-90,-82,-80,-78,-76,-74,-72,-70,-68,-66,-64,-62,-60,-58,-56,-54,-52,-48,\
    #               -46,-44,-42,-40,-38,-36,-34,-32,-30,-28,-26,-24,-22,-20,-18,-16,-14,-10]
    xmin_z = 8000; xmax_z = 10000; xint_z = 100# z interval
    color_vals = np.linspace(-1e-8, 1e-8, 41)

# color_cols = ['#ffffff','#ffe0ff','#ffc1ff','#ffa2ff','#ef7df4','#d854e6','#c22cd7','#a621d6',\
#               '#823ae4','#5f52f3','#416cfe','#3f8bf8','#3eaaf3','#3cc9ed','#51cdc0','#69cf8f',\
#               '#80d15d','#9dd947','#bee341','#deed3b','#f4ec36','#f8d533','#fcbe30','#ffa62c',\
#               '#ff8c1e','#ff7210','#ff5703','#ff3e36','#ff2579','#ff0bbb','#ff17df','#ff41eb',\
#               '#ff6af7','#ff92ff','#ffb7ff','#ffdbff','#ffffff']

vmin = np.min(np.min(color_vals))#levels)#-20#np.min(t)
vmax = np.max(np.max(color_vals)) #-5#np.max(t)

if region == 'Canada':
    shrinkscale = 0.3
else:
    shrinkscale = 0.6

cf = ax.contourf(LON, LAT, z1_del2, color_vals, cmap='RdBu_r', alpha=1, vmin=vmin, vmax=vmax)
cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Specific humidity [°C]', shrink=shrinkscale)
cbar.set_label("del2(z)", fontsize=16)  
cbar.ax.tick_params(labelsize=12) 

# ct = ax.contour(LON, LAT, z1_del2, levels=30, colors='white', linewidths=1)
# ax.clabel(ct, inline=True, fontsize=11, fmt='%d')#'%.2f')

cs = ax.contour(LON, LAT, z1, levels=np.arange(xmin_z,xmax_z+1,xint_z), colors='black', linewidths=1)
ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')

#Wind barbs are plotted every n gridpoint
# cuv = ax.barbs(LON[::n,::n], LAT[::n,::n], uu[::n,::n], vv[::n,::n])

if (pres == 300):
#     windmag = np.sqrt((uu/1.94384449)**2 + (vv/1.94384449)**2)
#     windmag = np.array(windmag)
#     windmag[windmag<30] = np.nan
#     windmag[windmag>=30] = 1
#     ax.contourf(LON, LAT, windmag, levels=[1,1.1], colors='none', alpha=0.2)
    plt.title('%d hPa Geopotential [gpdm, contour], del2(Z) (color)\n\
    %sUTC' %(pres,timt), fontsize=20)
else:
    plt.title('%d hPa Geopotential [gpdm, contour], del2(Z) (color)\n\
    %sUTC' %(pres,timt), fontsize=20)
if case == 'tornado':
    ax.scatter(lontstart, lattstart, s=150, marker="v", edgecolors='k', color='yellow')
    ax.scatter(lont, latt, s=150, marker="v", edgecolors='k', color='green', alpha=0.5)
else:
    ax.scatter(lontstart, lattstart, s=150, marker="^", edgecolors='k', color='cyan')
    ax.scatter(lont, latt, s=150, marker="^", edgecolors='k', color='green', alpha=0.5)

iyl,ixl = mc.find_peaks(z1, maxima=False)
iyh,ixh = mc.find_peaks(z1)
# for i in range(len(iyl)):
#     lonmin = lon1[ixl[i]]
#     latmin = lat1[iyl[i]]
#     ax.scatter(lonmin, latmin, s=400, marker='*', color='k', facecolor='b', linewidth=1.5)
# for j in range(len(iyh)):
#     lonmax = lon1[ixh[j]]
#     latmax = lat1[iyh[j]]
#     ax.scatter(lonmax, latmax, s=400, marker='*', color='k', facecolor='r', linewidth=1.5)

#%
il = list(zip(iyl,ixl))
ih = list(zip(iyh,ixh))
lll = list(zip(lat1[iyl], lon1[ixl]))
llh = list(zip(lat1[iyh], lon1[ixh]))

center_lat = lat1[80]
center_lon = lon1[140]
center = (center_lat, center_lon)

dist = np.zeros(shape=(len(lll),))
for k in range(len(lll)):
    dist[k] = np.sqrt(sum((a-b)**2 for a,b in zip(lll[k],center)))
    # print(f"lat/lon = {lll[k]} --> dist. from center = {dist[k]:.2f}")
ilow = np.argmin(dist)
latl = lll[ilow][0]
lonl = lll[ilow][1]
low = (latl, lonl)

# ax.scatter(lonl, latl, s=800, marker='*', color='k', facecolor='w', linewidth=1.5)
ax.text(lonl, latl, "L", color='b', fontsize=36, fontweight='bold', horizontalalignment='center', verticalalignment='center')

iy1,ix1 = mc.find_peaks(z1_del2)
lat2 = lat1[iy1][(lat1[iy1]<latl)]
lon2 = lon1[ix1][(lat1[iy1]<latl)]
lldel = list(zip(lat2,lon2))

# for i in range(len(lldel)):
#     ax.scatter(lldel[i][1], lldel[i][0], s=500, marker='.', color='k')

dist2 = np.zeros(shape=(len(lldel),))
for k in range(len(lldel)):
    dist2[k] = np.sqrt(sum((a-b)**2 for a,b in zip(lldel[k],low)))

idel = np.argmin(dist2)
latdel = lldel[idel][0]
londel = lldel[idel][1]

ax.scatter(londel, latdel, s=300, marker='o', color='k', facecolor='w', linewidth=1.5)

trough_angle = -1*np.arctan2(lonl-londel, latl-latdel)*180/np.pi

if (trough_angle >= -30) & (trough_angle <= 30):
    trough_tilt = "neutral"
elif (trough_angle > 30):
    trough_tilt = "negative"
elif (trough_angle < -30):
    trough_tilt = "positive"

print(f"{pres} mb trough axis = {trough_tilt}")

#%% Find trough tilt and loop through pressure levels

for pres in [500,300]:
    dat = data.sel(pressure_level=pres, valid_time=timt)
    z = dat["z"].values/9.81
    lon = dat["longitude"].values
    lat = dat["latitude"].values
    u = dat["u"].values
    v = dat["v"].values
    LON,LAT = np.meshgrid(lon,lat)


    iyl,ixl = mc.find_peaks(z, maxima=False)
    iyh,ixh = mc.find_peaks(z)

    dx,dy = mc.lat_lon_grid_deltas(lon,lat)

    # z_del2 = mc.geospatial_laplacian(dat["z"]/9.81, crs=ccrs.PlateCarree()).values
    # z_del2 = mc.smooth_gaussian(z_del2, 10)
    
    z2 = mc.smooth_gaussian(dat["z"]/9.81, 10)
    z_del2 = mc.geospatial_laplacian(z2, crs=ccrs.PlateCarree()).values
    


    #%
    fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()})
    ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle='-')
    ax.add_feature(cfeature.STATES, linestyle='-')
    ax.add_feature(cfeature.COASTLINE)

    if (pres == 850):
        xmin_z = 800; xmax_z = 1800; xint_z = 25     # z interval
        color_vals = np.linspace(-1e-8, 1e-8, 41)

    elif (pres == 700):
        xmin_z = 2400; xmax_z = 4000; xint_z = 50    # z interval
        color_vals = np.linspace(-1e-8, 1e-8, 41)

    elif (pres == 500):
        xmin_z = 4800; xmax_z = 6500; xint_z = 50    # z interval
        color_vals = np.linspace(-1e-8, 1e-8, 41)

    elif (pres == 300):
        xmin_z = 8000; xmax_z = 10000; xint_z = 100# z interval
        color_vals = np.linspace(-1e-8, 1e-8, 41)


    vmin = np.min(np.min(color_vals))#levels)#-20#np.min(t)
    vmax = np.max(np.max(color_vals)) #-5#np.max(t)

    if region == 'Canada':
        shrinkscale = 0.3
    else:
        shrinkscale = 0.6

    cf = ax.contourf(LON, LAT, z_del2, color_vals, cmap='RdBu_r', alpha=1, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Specific humidity [°C]', shrink=shrinkscale)
    cbar.set_label("del2(z)", fontsize=16)  
    cbar.ax.tick_params(labelsize=12) 

    # ct = ax.contour(LON, LAT, z_del2, levels=30, colors='white', linewidths=1)
    # ax.clabel(ct, inline=True, fontsize=11, fmt='%d')#'%.2f')

    cs = ax.contour(LON, LAT, z, levels=np.arange(xmin_z,xmax_z+1,xint_z), colors='black', linewidths=1)
    ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')

    #Wind barbs are plotted every n gridpoint
    # cuv = ax.barbs(LON[::n,::n], LAT[::n,::n], u[::n,::n]*1.94384449, v[::n,::n]*1.94384449)

    if (pres == 300):
    #     windmag = np.sqrt((u)**2 + (v)**2)
    #     windmag = np.array(windmag)
    #     windmag[windmag<30] = np.nan
    #     windmag[windmag>=30] = 1
    #     ax.contourf(LON, LAT, windmag, levels=[1,1.1], colors='none', alpha=0.2)
        plt.title('%d hPa Geopotential [gpdm, contour], del2(Z) (color)\n\
        %sUTC' %(pres,timt), fontsize=20)
    else:
        plt.title('%d hPa Geopotential [gpdm, contour], del2(Z) (color)\n\
        %sUTC' %(pres,timt), fontsize=20)
    if case == 'tornado':
        ax.scatter(lontstart, lattstart, s=150, marker="v", edgecolors='k', color='yellow')
        ax.scatter(lont, latt, s=150, marker="v", edgecolors='k', color='green', alpha=0.5)
    else:
        ax.scatter(lontstart, lattstart, s=150, marker="^", edgecolors='k', color='cyan')
        ax.scatter(lont, latt, s=150, marker="^", edgecolors='k', color='green', alpha=0.5)

    iyl,ixl = mc.find_peaks(z, maxima=False)
    iyh,ixh = mc.find_peaks(z)

    #%
    il = list(zip(iyl,ixl))
    ih = list(zip(iyh,ixh))
    lll = list(zip(lat[iyl], lon[ixl]))
    llh = list(zip(lat[iyh], lon[ixh]))

    center_lat = lat[80]
    center_lon = lon[140]
    center = (center_lat, center_lon)

    dist = np.zeros(shape=(len(lll),))
    for k in range(len(lll)):
        dist[k] = np.sqrt(sum((a-b)**2 for a,b in zip(lll[k],center)))
        # print(f"lat/lon = {lll[k]} --> dist. from center = {dist[k]:.2f}")
    ilow = np.argmin(dist)
    latl = lll[ilow][0]
    lonl = lll[ilow][1]
    low = (latl, lonl)

    # ax.scatter(lonl, latl, s=800, marker='*', color='k', facecolor='w', linewidth=1.5)
    ax.text(lonl, latl, "L", color='b', fontsize=36, fontweight='bold', horizontalalignment='center', verticalalignment='center')

    iy1,ix1 = mc.find_peaks(z_del2)
    lat2 = lat[iy1][(lat[iy1]<latl)]
    lon2 = lon[ix1][(lat[iy1]<latl)]
    lldel = list(zip(lat2,lon2))
    
    dist2 = np.zeros(shape=(len(lldel),))
    for k in range(len(lldel)):
        dist2[k] = np.sqrt(sum((a-b)**2 for a,b in zip(lldel[k],low)))

    idel = np.argmin(dist2)
    latdel = lldel[idel][0]
    londel = lldel[idel][1]

    ax.scatter(londel, latdel, s=300, marker='o', color='k', facecolor='w', linewidth=1.5)

    trough_angle = -1*np.arctan2(lonl-londel, latl-latdel)*180/np.pi

    if (trough_angle >= -30) & (trough_angle <= 30):
        trough_tilt = "neutral"
    elif (trough_angle > 30):
        trough_tilt = "negative"
    elif (trough_angle < -30):
        trough_tilt = "positive"

    print(f"{pres} mb trough axis = {trough_tilt}")
    
    












