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

fp = "C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/tor_outbreaks/"

case='tornado'
#case='hail'

# region='Alberta'
#region='SasMan'
region='GreatLakes'
# region='Canada'

yyyyt = 2026
mmt   = 6
ddt   = 30
hht   = 16

# yyyyt = 2021; mmt = 8; ddt = 11; tor hht = 20,21 (18-22)
# yyyyt = 2025; mmt = 6; ddt = 23-24; tor hht = 20,22,1 (19-23, 0-3)
# yyyyt = 2022; mmt = 5; ddt = 30-31; tor hht = 0,1,2 (22-23, 0-3)
# yyyyt = 2022; mmt = 5; ddt = 21; tor hht = 15,17 (13-18)
# yyyyt = 2025; mmt = 7; ddt = 24; bow/db hht = 21-23 (19-23)
# yyyyt = 2026; mmt = 4; ddt = 15; tor/null hht = 5,6 (2-8)
# yyyyt = 2026; mmt = 6; ddt = 30; tor hht = 16,17 (13-17)
# yyyyt = 2026; mmt = 7; ddt = 3; hht = 22-00 (20-23, 0-2)

timt="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))


# lattstart = 46.4797  #Latitude of tornado/beginning of hail swath (from NTP/NHP database)
# lontstart = -83.6403 #Longitude of tornado/beginning of hail swath (from NTP/NHP database)

if (yyyyt==2021) & (mmt==8) & (ddt==11):
    lattstart = [46.4797, 46.5449, 46.5724, 46.4657, 46.3870, 46.6557]
    lontstart = [-83.6403, -83.5907, -83.2889, -83.2155, -82.7493, -82.4617]
    figfolder = fp+"/figs/20210811/"
elif (yyyyt==2025) & (mmt==6) & ((ddt==23)|(ddt==24)):
    lattstart = [48.3975, 48.4647, 48.2033, 47.9149, 46.7872, 46.7757]
    lontstart = [-75.5880, -75.4900, -73.5234, -73.5056, -70.7827, -70.4072]
    figfolder = fp+"/figs/20250623/"
elif (yyyyt==2022) & (mmt==5) & ((ddt==30)|(ddt==31)):
    lattstart = [48.6148, 48.6699, 48.9114, 49.3689, 48.6757]
    lontstart = [-93.5304, -93.5219, -93.5778, -93.0754, -92.2308]
    figfolder = fp+"/figs/20220530/"
elif (yyyyt==2022) & (mmt==5) & (ddt==21):
    lattstart = [43.0179, 42.9217, 44.1058, 44.1755]
    lontstart = [-81.2216, -81.1977, -79.1458, -78.7722]
    figfolder = fp+"/figs/20220521/"
elif (yyyyt==2025) & (mmt==7) & (ddt==24):
    lattstart = [43.4700, 44.2500, 45.0000]
    lontstart = [-81.1800, -80.5000, -79.2500]
    figfolder = fp+"/figs/20250724/"
elif (yyyyt==2026) & (mmt==4) & (ddt==15):
    lattstart = [42.2718, 42.2684, 42.2393, 42.2423]
    lontstart = [-83.7524, -83.2104, -83.0576, -82.8484]
    figfolder = fp+"/figs/20260415/"
elif (yyyyt==2026) & (mmt==6) & (ddt==30):
    lattstart = [44.4640, 44.3021, 44.3303, 45.6315, 44.6890]
    lontstart = [-76.7289, -76.5386, -76.9268, -77.8432, -76.9664]
    figfolder = fp+"/figs/20260630/"
elif (yyyyt==2026) & (mmt==7) & ((ddt==3)|(ddt==4)):
    lattstart = [42.2393, 42.4226, 42.7505, 42.2766]
    lontstart = [-83.0576, -82.1339, -81.7004, -83.7341]
    figfolder = fp+"/figs/20260703/"


# latt = 46.50       # latitude of sounding
# lont = -83.75     # longitude of sounding #if you do not need this, uncomment the next two lines




# figfolder = fp+f"figs/"
data_preslevs = fp+f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_preslevs.nc"
data_singlevs = fp+f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_singlevs.nc"



#% Read data files

# Pressure level data
data = xr.open_dataset(data_preslevs)
# data.head(5)

# Single level data
datas = xr.open_dataset(data_singlevs)
# datas.head(5)


#%% Make plots

from era5utils import *

#% Define plot domain - can add other regions here

#Plotdomain:
#lonW=-125
#lonE=-90
#latS=45
#latN=65

# region='Alberta'
#region='SasMan'
# region='GreatLakes'
region='Canada'

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



data00 = data.sel(valid_time=timt)
z00 = data00['z'].values/9.81
u00 = data00['u'].values*units('m/s')
v00 = data00['v'].values*units('m/s')
p00 = data00.pressure_level.values * units.hPa
lat = data00['latitude'].values
lon = data00['longitude'].values

lati = np.argmin(np.abs(lat-np.mean(lattstart)))
loni = np.argmin(np.abs(lon-np.mean(lontstart)))
latt = lat[lati]
lont = lon[loni]




figfolder=None

PlotCAPEMap(data,datas,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,pngfolder=figfolder)

PlotPressureMaps(data,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,pngfolder=figfolder)

PlotMoistureMap2(data,datas,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,pngfolder=figfolder)

PlotSRHMaps(data,datas,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,pngfolder=figfolder)


# latt=51.5  #Latitude of tornado/beginning of hail swath (from NTP/NHP database)
# lont=-114.0

# lati=np.where(lat==latt)
# loni=np.where(lon==lont)



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
    
    

#%% Plotting my own upper air/surface/sounding parameter maps

fp = "C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/tor_outbreaks/"

# region = 'Alberta'
# region = 'SasMan'
# region = 'GreatLakes'
region = 'Canada'
# region = None

yyyyt = 2026
mmt   = 7
ddt   = 3
hht   = 22



figsave = False

# yyyyt = 2021; mmt = 8; ddt = 11; tor hht = 20,21 (17-22)
# yyyyt = 2025; mmt = 6; ddt = 23-24; tor hht = 20,22,1 (17-23, 0-3)
# yyyyt = 2022; mmt = 5; ddt = 30-31; tor hht = 0,1,2 (21-23, 0-3)
# yyyyt = 2022; mmt = 5; ddt = 21; tor hht = 15,17 (12-18)
# yyyyt = 2026; mmt = 6; ddt = 30; tor hht = 16,17 (13-17)
# yyyyt = 2025; mmt = 7; ddt = 24-25; bow/db hht = 21-23 (19-23, 0-1)
# yyyyt = 2026; mmt = 7; ddt = 3-4; hht = 22-00 (19-23, 0-2)
# yyyyt = 2026; mmt = 4; ddt = 15; tor/null hht = 5,6 (2-8)

timt="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))


# lattstart = 46.4797  #Latitude of tornado/beginning of hail swath (from NTP/NHP database)
# lontstart = -83.6403 #Longitude of tornado/beginning of hail swath (from NTP/NHP database)

if (yyyyt==2021) & (mmt==8) & (ddt==11):
    lattstart = [46.4797, 46.5449, 46.5724, 46.4657, 46.3870, 46.6557]
    lontstart = [-83.6403, -83.5907, -83.2889, -83.2155, -82.7493, -82.4617]
    figfolder = fp+"/figs/20210811/"
elif (yyyyt==2025) & (mmt==6) & ((ddt==23)|(ddt==24)):
    lattstart = [48.3975, 48.4647, 48.2033, 47.9149, 46.7872, 46.7757]
    lontstart = [-75.5880, -75.4900, -73.5234, -73.5056, -70.7827, -70.4072]
    figfolder = fp+"/figs/20250623/"
elif (yyyyt==2022) & (mmt==5) & ((ddt==30)|(ddt==31)):
    lattstart = [48.6148, 48.6699, 48.9114, 49.3689, 48.6757]
    lontstart = [-93.5304, -93.5219, -93.5778, -93.0754, -92.2308]
    figfolder = fp+"/figs/20220530/"
elif (yyyyt==2022) & (mmt==5) & (ddt==21):
    lattstart = [43.0179, 42.9217, 44.1058, 44.1755]
    lontstart = [-81.2216, -81.1977, -79.1458, -78.7722]
    figfolder = fp+"/figs/20220521/"
elif (yyyyt==2025) & (mmt==7) & (ddt==24):
    lattstart = [43.4700, 44.2500, 45.0000]
    lontstart = [-81.1800, -80.5000, -79.2500]
    figfolder = fp+"/figs/20250724/"
elif (yyyyt==2026) & (mmt==4) & (ddt==15):
    lattstart = [42.2718, 42.2684, 42.2393, 42.2423]
    lontstart = [-83.7524, -83.2104, -83.0576, -82.8484]
    figfolder = fp+"/figs/20260415/"
elif (yyyyt==2026) & (mmt==6) & (ddt==30):
    lattstart = [44.4640, 44.3021, 44.3303, 45.6315, 44.6890]
    lontstart = [-76.7289, -76.5386, -76.9268, -77.8432, -76.9664]
    figfolder = fp+"/figs/20260630/"
elif (yyyyt==2026) & (mmt==7) & ((ddt==3)|(ddt==4)):
    lattstart = [42.2393, 42.4226, 42.7505, 42.2766]
    lontstart = [-83.0576, -82.1339, -81.7004, -83.7341]
    figfolder = fp+"/figs/20260703/"



[lonW,lonE,latS,latN] = [np.median(lontstart)-13, np.median(lontstart)+13, np.median(lattstart)-7, np.median(lattstart)+7]

if region=='GreatLakes':    
    [lonW,lonE,latS,latN] = [-93.0,-70.0,39.5,52.0]
    n=5
    shrinkscale = 0.6
elif region=='SasMan':
    [lonW,lonE,latS,latN] = [-114.0,-88.0,46.5,60.5]
    n=5
    shrinkscale = 0.6
elif region=='Alberta':      
    [lonW,lonE,latS,latN] = [-128.0,-102.0,46.5,60.5]
    n=5
    shrinkscale = 0.6
elif region=='Canada':
    [lonW,lonE,latS,latN] = [-132.0,-58.0,35.0,75.0]
    n=10
    shrinkscale = 0.3

# figfolder = fp+'figs/'
data_preslevs = fp+f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_preslevs.nc"
data_singlevs = fp+f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_singlevs.nc"

data = xr.open_dataset(data_preslevs)
datas = xr.open_dataset(data_singlevs)

# def PlotMap(data, datas, plotvar, pressure_level, timt, lontstart, lattstart, region="Canada", lonlatWESN=None):
### plotvar ###
### On pressure levels-
# z:  Geopotential (m2/s2)
# pv: Potential vorticity (K m2/kg/s)
# r:  Relative humidity (%)
# q:  Specific humidity (kg/kg)
# t:  Temperature (K)
# u:  U wind (m/s)
# v:  V wind (m/s)
# w:  Vertical velocity (Pa/s)
# vo: Relative vorticity (1/s)
### On single levels-
# u10:  10-m U wind (m/s)
# v10:  10-m V wind (m/s)
# d2m:  2-m dewpoint (K)
# t2m:  2-m temperature (K)
# sp:   Surface pressure (Pa)
# u100: 100-m U wind (m/s)
# v100: 100-m V wind (m/s)
# cape: CAPE (J/kg)
# cin:  CIN (J/kg)
# z:    Geopotential of surface (m2/s2)
# tcw:  Total column water (kg/m2)


dats = datas.sel(valid_time=timt)
orog = dats['z'].values/9.81 #Geopotential height of surface (m)
u10 = dats['u10'].values #10-m u wind (m/s)
v10 = dats['v10'].values #10-m v wind (m/s)
uu10 = u10*1.94384449 #10-m u wind to kts
vv10 = v10*1.94384449 #10-m v wind to kts
u100 = dats['u100'].values #100-m u wind (m/s)
v100 = dats['v100'].values #100-m v wind (m/s)
uu100 = u100*1.94384449 #100-m u wind to kts
vv100 = v100*1.94384449 #100-m v wind to kts
d2m = dats['d2m'].values-273.15 #2-m dewpoint (C)
t2m = dats['t2m'].values-273.15 #2-m temperature (C)
psfc = dats['sp'].values/100 #Surface pressure (hPa)
cape = dats['cape'].values #CAPE
cin = dats['cin'].values #CIN
# mslp = psfc * (1 - (0.0065*orog)/(t2m+273.15))**(-1*(9.81*0.0289644)/(8.31447*0.0065)) # from https://absolutepressurecalculator.com/how-to-calculate-mean-sea-level-pressure.php
mslp = dats['msl'].values/100 #MSLP (hPa)


datt = data.sel(valid_time=timt)
prs = datt.pressure_level.values * units.hPa
lat = datt['latitude'].values
lon = datt['longitude'].values
X,Y = np.meshgrid(lon,lat)

z = datt['z'].values/9.81 #Geopotential height (m)
t = datt['t'].values-273.15 #Temperature (C)
rh = datt['r'].values*units('percent') #Relative humidity
q = datt['q'].values*units('kg/kg') #Specific humidity
u = datt['u'].values #U wind (m/s)
v = datt['v'].values #V wind (m/s)
omega = datt['w'].values*units('Pa/s') #omega (Pa/s)
w = mc.vertical_velocity(omega, np.moveaxis(np.tile(prs, (len(lat),len(lon),1)), -1, 0), t*units.degC, q)
pv = datt['pv'].values
vort = datt['vo'].values

lati = np.argmin(np.abs(lat-np.mean(lattstart)))
loni = np.argmin(np.abs(lon-np.mean(lontstart)))
latt = lat[lati]
lont = lon[loni]


data.close()
datas.close()
datt.close()
dats.close()

#%% Upper air maps

colors = ['dodgerblue','lightskyblue','cyan','mediumpurple','blueviolet','mediumvioletred']
spc_wspd = ListedColormap(colors, name="spc_wspd")

# plot_params = {
#     'z850':    {'cm': 'LangRainbow12',  'label': "850 mb geopotential hgt (m)",  'levels':np.arange(900,1830,30)},
#     'z700':    {'cm': 'LangRainbow12',  'label': "700 mb geopotential hgt (m)",  'levels':np.arange(2100,3630,30)},
#     'z500':    {'cm': 'LangRainbow12',  'label': "500 mb geopotential hgt (m)",  'levels':np.arange(4500,6630,60)},
#     'z300':    {'cm': 'LangRainbow12',  'label': "300 mb geopotential hgt (m)",  'levels':np.arange(7500,9930,60)},
#     't850':    {'cm': 'LangRainbow12',  'label': "850 mb temperature (C)",   'levels':np.arange(-20,42,2)},
#     't700':    {'cm': 'LangRainbow12',  'label': "850 mb temperature (C)",   'levels':np.arange(-30,32,2)},
#     't500':    {'cm': 'LangRainbow12',  'label': "850 mb temperature (C)",   'levels':np.arange(-40,22,2)},
#     'wspd':    {'cm':  spc_wspd,        'label': "Wind speed (m s$^{-1}$)",  'levels':np.arange(40,180,20)},
#     'temp':    {'cm': 'HomeyerRainbow', 'label': "Temperature (C)"},
#     'dewpt':   {'cm': 'YlGnBu',         'label': "Dewpoint (C)"},
#     'rh':      {'cm': 'YlGn',           'label': "Relative humidity (%)"},
#     'q':       {'cm': 'YlGn',           'label': "q (g kg$^{-1}$)"},
#     'vort':    {'cm': 'HomeyerRainbow', 'label': "\u03B6 (s$^{-1}$)"},
#     'cape':    {'cm': 'LangRainbow12',  'label': "CAPE (J kg$^{-1}$)"},
#     'srh':     {'cm': 'LangRainbow12',  'label': "SRH (m$^{2}$ s$^{-2}$)"}
# }


# figsave = True




### Upper air maps ###
for presl in [300,500,700,850]:
    
    zz = datt.sel(pressure_level=presl)['z'].values/9.81
    uu = datt.sel(pressure_level=presl)['u'].values*1.94384449
    vv = datt.sel(pressure_level=presl)['v'].values*1.94384449
    wspd = np.sqrt(uu**2 + vv**2)
    tt = datt.sel(pressure_level=presl)['t'].values-273.15
    qq = datt.sel(pressure_level=presl)['q'].values
    e = (qq * presl) / (0.622+qq)
    td = 243.5 / ((17.67/(np.log(e/6.112)))-1)
    
    
    fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()}, layout='tight')
    ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle='-')
    ax.add_feature(cfeature.STATES, linestyle='-')
    ax.add_feature(cfeature.COASTLINE)
    
    # Geopotential height, wind speed and wind barbs (kts)
    if presl == 300:
        cf = ax.contourf(X, Y, wspd, levels=np.linspace(60,180,7), cmap=spc_wspd, vmin=60, vmax=180, alpha=0.35)
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Wind speed [kts]', shrink=shrinkscale)
        cbar.set_label('Wind speed [kts]', fontsize=16)  
        cbar.ax.tick_params(labelsize=12)
        cw = ax.contour(X, Y, wspd, levels=np.linspace(60,160,6), colors='b', linestyles='-', linewidths=1)
        ax.clabel(cw, inline=True, fontsize=11, fmt='%d')#.1f')
        
        cs = ax.contour(X, Y, zz, levels=np.arange(7500,9930,60), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
        
        ax.set_title(f"300-mb wind speed (shaded), geopotential height\n {timt[:10]} {timt[11:13]}Z ", fontsize=20)
    
    # Geopotential height, wind speed and wind barbs (kts)
    if presl == 500:
        cf = ax.contourf(X, Y, wspd, levels=np.linspace(40,160,7), colors=colors, vmin=40, vmax=160, alpha=0.35)
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Wind speed [kts]', shrink=shrinkscale)
        cbar.set_label('Wind speed [kts]', fontsize=16)  
        cbar.ax.tick_params(labelsize=12)
        cw = ax.contour(X, Y, wspd, levels=np.linspace(40,140,6), colors='b', linestyles='-', linewidths=1)
        ax.clabel(cw, inline=True, fontsize=11, fmt='%d')#.1f')
        
        cs = ax.contour(X, Y, zz, levels=np.arange(4500,6630,60), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
        
        ax.set_title(f"500-mb wind speed (shaded), geopotential height\n {timt[:10]} {timt[11:13]}Z ", fontsize=20)
    
    # Geopotential height, relative humidity, 700-500 mb lapse rates, wind barbs (kts)
    if presl == 700:
        rr = np.mean(datt.sel(pressure_level=slice(700,500))['r'].values, axis=0)
        t75 = datt.sel(pressure_level=slice(700,500))['t'].values-273.15
        z75 = datt.sel(pressure_level=slice(700,500))['z'].values/9.81/1000
        lr = -1*(t75[-1,:,:]-t75[0,:,:])/(z75[-1,:,:]-z75[0,:,:])
        
        cf = ax.contourf(X, Y, rr, levels=[70,101], colors=['mediumaquamarine'], alpha=0.4)
        # cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Rel. humidity [%]', shrink=shrinkscale)
        # cbar.set_label('Rel. humidity [%]', fontsize=16)  
        # cbar.ax.tick_params(labelsize=12)
        cr = ax.contour(X, Y, rr, levels=np.arange(70,101,10), colors='teal', linestyles='-', linewidths=1)
        ax.clabel(cr, inline=True, fontsize=11, fmt='%d')#.1f')
        # ct1 = ax.contour(X, Y, tt, levels=np.arange(2,10,2), colors='r', linestyles='--', linewidths=1.5)
        # ct2 = ax.contour(X, Y, tt, levels=np.arange(10,30,2), colors='r', linestyles='--', linewidths=2.5)
        # ct3 = ax.contour(X, Y, tt, levels=np.arange(-30,2,2), colors='b', linestyles='--', linewidths=1.5)
        # ax.clabel(ct1, inline=True, fontsize=11, fmt='%d')#.1f')
        # ax.clabel(ct2, inline=True, fontsize=11, fmt='%d')#.1f')
        # ax.clabel(ct3, inline=True, fontsize=11, fmt='%d')#.1f')
        cf2 = ax.contourf(X, Y, lr, levels=[8,10], colors=['r'], alpha=0.3)
        cl = ax.contour(X, Y, lr, levels=np.linspace(6,9,7), colors='r', linewidths=[1,1.5,1,1,1,1,2.5], linestyles=['--','--','-','-','-','-','-'])
        ax.clabel(cl, inline=True, fontsize=11, fmt='%.1f')#.1f')
        
        cs = ax.contour(X, Y, zz, levels=np.arange(2100,3630,30), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
        
        ax.set_title(f"700-mb RH>70% (green,shaded), geopotential height (black), winds (barbs), 700-500-mb lapse rate (red)\n {timt[:10]} {timt[11:13]}Z ", fontsize=20)
    
    # Geopotential height, temperature, dewpoint, wind barbs (kts)
    if presl == 850:
        if region == 'Canada':
            dlevs = np.arange(10,30,5); tlevs1 = np.arange(5,20,5); tlevs2 = np.arange(20,40,5); tlevs3 = np.arange(-30,5,5)
        else:
            dlevs = np.arange(8,30,2); tlevs1 = np.arange(2,20,2); tlevs2 = np.arange(20,40,2); tlevs3 = np.arange(-30,2,2)
        
        cf = ax.contourf(X, Y, td, levels=dlevs, colors=['mediumaquamarine'], alpha=0.4)
        # cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Dewpoint [C]', shrink=shrinkscale)
        # cbar.set_label('Dewpoint [C]', fontsize=16)  
        # cbar.ax.tick_params(labelsize=12)
        cd = ax.contour(X, Y, td, levels=dlevs, colors='teal', linestyles='-', linewidths=1)
        ct1 = ax.contour(X, Y, tt, levels=tlevs1, colors='r', linestyles='--', linewidths=1.5)
        ct2 = ax.contour(X, Y, tt, levels=tlevs2, colors='r', linestyles='--', linewidths=2)
        ct3 = ax.contour(X, Y, tt, levels=tlevs3, colors='b', linestyles='--', linewidths=1.5)
        ax.clabel(cd, inline=True, fontsize=11, fmt='%d')#.1f')
        ax.clabel(ct1, inline=True, fontsize=11, fmt='%d')#.1f')
        ax.clabel(ct2, inline=True, fontsize=11, fmt='%d')#.1f')
        ax.clabel(ct3, inline=True, fontsize=11, fmt='%d')#.1f')
        
        cs = ax.contour(X, Y, zz, levels=np.arange(900,1830,30), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
        
        ax.set_title(f"850-mb temperature (red/blue), dewpoint (green,shaded), geopotential height (black), winds (barbs)\n {timt[:10]} {timt[11:13]}Z ", fontsize=20)
    
    #Wind barbs are plotted every n gridpoint
    cuv = ax.barbs(X[::n,::n], Y[::n,::n], uu[::n,::n], vv[::n,::n])
    
    s = ax.scatter(lontstart, lattstart, s=250, marker="v", edgecolors='k', color='yellow', linewidth=2, zorder=10)
    ax.legend(handles=[s], labels=['Tornado'], loc='upper right', fontsize=20, bbox_to_anchor=(1,1.07), frameon=False)
    
    
    if region == 'Canada':
        figname = f"upperair_{presl}mb_{hht:02.0f}Z_Canada.png"
    else:
        figname = f"upperair_{presl}mb_{hht:02.0f}Z.png"
    
    if figsave:
        plt.savefig(figfolder+figname, dpi=300)
    # plt.show()



#%% Surface maps


# 2-m temperature, surface pressure, 10-m wind barbs (kts)
fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()}, layout='tight')
ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linestyle='-')
ax.add_feature(cfeature.STATES, linestyle='-')
ax.add_feature(cfeature.COASTLINE)

cf = ax.contourf(X, Y, t2m, levels=np.arange(-30,42,2), cmap='LangRainbow12')
cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Temperature [C]', shrink=shrinkscale)
cbar.set_label('Temperature [C]', fontsize=16)  
cbar.ax.tick_params(labelsize=12)
cp = ax.contour(X, Y, mslp, levels=np.arange(900,1040,4), colors='k', linewidths=2, linestyles='-')
ax.clabel(cp, inline=True, fontsize=11, fmt='%d')#.1f')
# cuv100 = ax.barbs(X[::n,::n], Y[::n,::n], uu100[::n,::n], vv100[::n,::n], color='r')
cuv = ax.barbs(X[::n,::n], Y[::n,::n], uu10[::n,::n], vv10[::n,::n], color='k')

s = ax.scatter(lontstart, lattstart, s=250, marker="v", edgecolors='k', color='yellow', linewidth=2, zorder=10)
ax.legend(handles=[s], labels=['Tornado'], loc='upper right', fontsize=20, bbox_to_anchor=(1,1.07), frameon=False)

ax.set_title(f"2-m temperature, MSLP (black), 10-m winds (barbs)\n {timt[:10]} {timt[11:13]}Z ", fontsize=20)

if region == 'Canada':
    figname = f"sfc_temp_{hht:02.0f}Z_Canada.png"
else:
    figname = f"sfc_temp_{hht:02.0f}Z.png"
    
if figsave:
    plt.savefig(figfolder+figname, dpi=300)
# plt.show()


# 2-m dewpoint, surface pressure, 10-m wind barbs
fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()}, layout='tight')
ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linestyle='-')
ax.add_feature(cfeature.STATES, linestyle='-')
ax.add_feature(cfeature.COASTLINE)

cf = ax.contourf(X, Y, d2m, levels=np.arange(-10,32,2), cmap='terrain_r')
cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Dewpoint [C]', shrink=shrinkscale)
cbar.set_label('Dewpoint [C]', fontsize=16)  
cbar.ax.tick_params(labelsize=12)
cp = ax.contour(X, Y, mslp, levels=np.arange(900,1040,4), colors='k', linewidths=2, linestyles='-')
ax.clabel(cp, inline=True, fontsize=11, fmt='%d')#.1f')
# cuv100 = ax.barbs(X[::n,::n], Y[::n,::n], uu100[::n,::n], vv100[::n,::n], color='r')
cuv = ax.barbs(X[::n,::n], Y[::n,::n], uu10[::n,::n], vv10[::n,::n], color='k')

s = ax.scatter(lontstart, lattstart, s=250, marker="v", edgecolors='k', color='yellow', linewidth=2, zorder=10)
ax.legend(handles=[s], labels=['Tornado'], loc='upper right', fontsize=20, bbox_to_anchor=(1,1.07), frameon=False)

ax.set_title(f"2-m dewpoint, MSLP (black), 10-m winds (barbs)\n {timt[:10]} {timt[11:13]}Z ", fontsize=20)

if region == 'Canada':
    figname = f"sfc_dewpt_{hht:02.0f}Z_Canada.png"
else:
    figname = f"sfc_dewpt_{hht:02.0f}Z.png"

if figsave:
    plt.savefig(figfolder+figname, dpi=300)
# plt.show()



# MSLP, 2-m temperature, 2-m dewpoint, 10-m wind barbs
fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()}, layout='tight')
ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linestyle='-')
ax.add_feature(cfeature.STATES, linestyle='-')
ax.add_feature(cfeature.COASTLINE)

# cf = ax.contourf(X, Y, d2m-t2m, levels=[-2,1], colors=['limegreen'], alpha=0.2)
ct1 = ax.contour(X, Y, t2m, levels=np.arange(-30,5,5), colors='b', linewidths=1.5, linestyles='--')
ax.clabel(ct1, inline=True, fontsize=11, fmt='%d')#.1f')
ct2 = ax.contour(X, Y, t2m, levels=np.arange(5,20,5), colors='r', linewidths=1.5, linestyles='--')
ax.clabel(ct2, inline=True, fontsize=11, fmt='%d')#.1f')
ct3 = ax.contour(X, Y, t2m, levels=np.arange(20,45,5), colors='r', linewidths=1.5, linestyles='-')
ax.clabel(ct3, inline=True, fontsize=11, fmt='%d')#.1f')
cd1 = ax.contour(X, Y, d2m, levels=np.arange(5,20,5), colors='teal', linewidths=1.5, linestyles='--')
ax.clabel(cd1, inline=True, fontsize=11, fmt='%d')#.1f')
cd2 = ax.contour(X, Y, d2m, levels=np.arange(20,45,5), colors='teal', linewidths=1.5, linestyles='-')
ax.clabel(cd2, inline=True, fontsize=11, fmt='%d')#.1f')
cp = ax.contour(X, Y, mslp, levels=np.arange(900,1040,4), colors='k', linewidths=2, linestyles='-')
ax.clabel(cp, inline=True, fontsize=11, fmt='%d')#.1f')
cuv = ax.barbs(X[::n,::n], Y[::n,::n], uu10[::n,::n], vv10[::n,::n], color='k')

s = ax.scatter(lontstart, lattstart, s=250, marker="v", edgecolors='k', color='yellow', linewidth=2, zorder=10)
ax.legend(handles=[s], labels=['Tornado'], loc='upper right', fontsize=20, bbox_to_anchor=(1,1.07), frameon=False)

ax.set_title(f"MSLP (black), 2-m temperature (red/blue), 2-m dewpoint (green), 10-m winds (barbs)\n {timt[:10]} {timt[11:13]}Z ", fontsize=20)

if region == 'Canada':
    figname = f"sfc_pres_{hht:02.0f}Z_Canada.png"
else:
    figname = f"sfc_pres_{hht:02.0f}Z.png"

if figsave:
    plt.savefig(figfolder+figname, dpi=300)
# plt.show()




#%% Composite parameter maps

zz = datt.sel(pressure_level=500)['z'].values/9.81
uu = u*1.94384449 #U wind to kts
vv = v*1.94384449 #V wind to kts
uu6000 = InterpolateToHeightAboveGround(z, orog, uu, 6000)
vv6000 = InterpolateToHeightAboveGround(z, orog, vv, 6000)

rr = np.mean(datt.sel(pressure_level=slice(700,500))['r'].values, axis=0)
t75 = datt.sel(pressure_level=slice(700,500))['t'].values-273.15
z75 = datt.sel(pressure_level=slice(700,500))['z'].values/9.81/1000
lr = -1*(t75[-1,:,:]-t75[0,:,:])/(z75[-1,:,:]-z75[0,:,:])


# CAPE, CIN>100 J/kg, 500 mb gph, 10-6000 m wind difference (barbs, kts)
fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()}, layout='tight')
ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linestyle='-')
ax.add_feature(cfeature.STATES, linestyle='-')
ax.add_feature(cfeature.COASTLINE)

colormap = cmocean.cm.thermal.copy()
colormap.set_over("gray")
colormap.set_under(alpha=0)
vmin = 100
vmax = 3000
levels = np.arange(vmin,vmax+1,100)
nlevels = len(levels)

norm = BoundaryNorm(np.linspace(vmin, vmax, nlevels+1), ncolors=colormap.N)

# cf = ax.contourf(X, Y, cape, levels=np.arange(vmin,vmax+1,100), cmap=colormap)
cf = ax.contourf(X, Y, cape, levels=np.arange(vmin,vmax+1,100), cmap=colormap, norm=norm, vmin=vmin, vmax=vmax, extend='both')
cbar = plt.colorbar(cf, ax=ax, orientation='vertical', shrink=shrinkscale)
cbar.set_label('CAPE [J/kg]', fontsize=16)  
cbar.ax.tick_params(labelsize=12)
ax.contourf(X, Y, cin, levels=[100,1000], colors=['w'], alpha=0.4)
# ax.contourf(X, Y, cin, levels=[25,100], colors=['w'], alpha=0.2)
cc = ax.contour(X, Y, cin, levels=[100], colors='w', linewidths=1, linestyles='-')
# ax.clabel(cc, inline=True, fontsize=11, fmt='%d')
cf3 = ax.contourf(X, Y, lr, levels=[8,10], colors=['r'], alpha=0.3)
cl = ax.contour(X, Y, lr, levels=[6,6.5,7,7.5,8,9], colors='r', linewidths=[1,1.5,1,1.5,1.5,2], linestyles=['--','--','-','-','-','-'])
ax.clabel(cl, inline=True, fontsize=10, fmt='%.1f')#.1f')

cs = ax.contour(X, Y, zz, levels=np.arange(4500,6630,60), colors='k', linewidths=2)#'white', linewidths=1)
ax.clabel(cs, inline=True, fontsize=11, fmt='%d')

cuv = ax.barbs(X[::n,::n], Y[::n,::n], uu6000[::n,::n]-uu10[::n,::n], vv6000[::n,::n]-vv10[::n,::n], color="k")

s = ax.scatter(lontstart, lattstart, s=250, marker="v", edgecolors='k', color='yellow', linewidth=2, zorder=10)
ax.legend(handles=[s], labels=['Tornado'], loc='upper right', fontsize=20, bbox_to_anchor=(1,1.07), frameon=False)

ax.set_title(f"CAPE, CIN (white), 0-6 km shear (barbs), 500-mb geopotential height (black), 700-500 mb lapse (red)\n {timt[:10]} {timt[11:13]}Z ", fontsize=20)

if region == 'Canada':
    figname = f"cape_cin_{hht:02.0f}Z_Canada.png"
else:
    figname = f"cape_cin_{hht:02.0f}Z.png"

if figsave:
    plt.savefig(figfolder+figname, dpi=300)
# plt.show()



#% 0-1 and 0-3 SRH, 500 mb gph, CAPE>200 and 1000 J/kg, 0-1 and 0-3 wind difference (barbs, kts)
for zdep in [1000,3000]:
    u2 = InterpolateToHeightAboveGround(z, orog, u, zdep)
    v2 = InterpolateToHeightAboveGround(z, orog, v, zdep)
    uu2 = u2*1.94384449
    vv2 = v2*1.94384449
    
    [srhpos,srhneg,srhtot] = getSRH(z, orog, u, v, zdep, u10, v10, storm_vector='rightmover')
    srh = srhtot.magnitude
    
    
    fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()}, layout='tight')
    ax.set_extent([lonW,lonE,latS,latN], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle='-')
    ax.add_feature(cfeature.STATES, linestyle='-')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAKES, alpha=0.2)
    
    if zdep == 1000:
        vmin = 10    #np.min(srh)
        vmax = 300   #np.max(srh)
        # nlevels = 30
    elif zdep == 3000:
        vmin = 10    #np.min(srh)
        vmax = 300   #np.max(srh)
        # nlevels = 30
    
    colormap = pyart.graph.cmweather.cm.LangRainbow12.copy()
    colormap.set_over("gray")
    colormap.set_under(alpha=0)

    if region == 'Canada':
        shrinkscale = 0.3
        scale = 5
    else:
        shrinkscale = 0.6
        scale = 10

    levels = np.arange(vmin,vmax+1,10)
    nlevels = len(levels)
    
    norm = BoundaryNorm(np.linspace(vmin, vmax, nlevels+1), ncolors=colormap.N)
    
    cf = ax.contourf(X, Y, srh, levels=levels, cmap=colormap, norm=norm, vmin=vmin, vmax=vmax, extend='both')
    cbar = plt.colorbar(cf, ax=ax, orientation='vertical', shrink=shrinkscale)
    cbar.set_label('SRH [m2/s2]', fontsize=16)  
    cbar.ax.tick_params(labelsize=12)
    
    cs = ax.contour(X, Y, zz, levels=np.arange(4500,6630,60), colors='k', linewidths=2)#'white', linewidths=1)
    ax.clabel(cs, inline=True, fontsize=11, fmt='%d')
    
    cc = ax.contour(X, Y, cape, levels=[200,1000], colors='b', linewidths=1.5, linestyles=['--','-'])
    ax.clabel(cc, inline=True, fontsize=11, fmt='%d')
    
    #Wind barbs are plotted every n gridpoint
    # cuv = ax.quiver(X[::n,::n], Y[::n,::n], u2[::n,::n], v2[::n,::n], angles='xy', scale_units='xy', scale=scale)
    # ax.quiverkey(cuv, X=1.08, Y=1.02, U=10, label="10 m/s", labelpos='N', fontproperties=dict(size=12))
    
    cuv = ax.barbs(X[::n,::n], Y[::n,::n], uu2[::n,::n]-uu10[::n,::n], vv2[::n,::n]-vv10[::n,::n], color="k")
    
    s = ax.scatter(lontstart, lattstart, s=250, marker="v", edgecolors='k', color='yellow', linewidth=2, zorder=10)
    ax.legend(handles=[s], labels=['Tornado'], loc='upper right', fontsize=20, bbox_to_anchor=(1,1.07), frameon=False)
    
    ax.set_title(f"0-{zdep/1000:.0f} km SRH, CAPE>200,1000 J/kg (blue), 0-{zdep/1000:.0f} km shear (barbs), 500-mb geopotential height (black)\n {timt[:10]} {timt[11:13]}Z ", fontsize=20)
    
    if region == 'Canada':
        figname = f"srh_0-{zdep/1000:.0f}km_{hht:02.0f}Z_Canada.png"
    else:
        figname = f"srh_0-{zdep/1000:.0f}km_{hht:02.0f}Z.png"
    
    if figsave:
        plt.savefig(figfolder+figname, dpi=300)
    # plt.show()



#%% 3-h height falls, vorticity

fp = "C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/tor_outbreaks/"

# region = 'Alberta'
# region = 'SasMan'
# region = 'GreatLakes'
region = 'Canada'
# region = None

yyyy = 2026
mm   = 4
dd   = 15
# hht   = 8


dbfile = open(fp+'tornado_locs.pkl', 'rb')
locs = pickle.load(dbfile)
tor1 = locs[f"{yyyy}{mm:02.0f}{dd:02.0f}"]['loc1']
dbfile.close

timestr = tor1['time_utc']
[yyyyt,mmt,ddt] = tor1['date_ymd']
hht = int(timestr[:2])


figsave = False

# yyyyt = 2021; mmt = 8; ddt = 11; tor hht = 20,21 (17-22)
# yyyyt = 2025; mmt = 6; ddt = 23-24; tor hht = 20,22,1 (17-23, 0-3)
# yyyyt = 2022; mmt = 5; ddt = 30-31; tor hht = 0,1,2 (21-23, 0-3)
# yyyyt = 2022; mmt = 5; ddt = 21; tor hht = 15,17 (12-18)
# yyyyt = 2026; mmt = 6; ddt = 30; tor hht = 16,17 (13-17)
# yyyyt = 2025; mmt = 7; ddt = 24-25; bow/db hht = 21-23 (19-23, 0-1)
# yyyyt = 2026; mmt = 7; ddt = 3-4; hht = 22-00 (19-23, 0-2)
# yyyyt = 2026; mmt = 4; ddt = 15; tor/null hht = 5,6 (2-8)




# lattstart = 46.4797  #Latitude of tornado/beginning of hail swath (from NTP/NHP database)
# lontstart = -83.6403 #Longitude of tornado/beginning of hail swath (from NTP/NHP database)

if (yyyyt==2021) & (mmt==8) & (ddt==11):
    lattstart = [46.4797, 46.5449, 46.5724, 46.4657, 46.3870, 46.6557]
    lontstart = [-83.6403, -83.5907, -83.2889, -83.2155, -82.7493, -82.4617]
    figfolder = fp+"/figs/20210811/"
elif (yyyyt==2025) & (mmt==6) & ((ddt==23)|(ddt==24)):
    lattstart = [48.3975, 48.4647, 48.2033, 47.9149, 46.7872, 46.7757]
    lontstart = [-75.5880, -75.4900, -73.5234, -73.5056, -70.7827, -70.4072]
    figfolder = fp+"/figs/20250623/"
elif (yyyyt==2022) & (mmt==5) & ((ddt==30)|(ddt==31)):
    lattstart = [48.6148, 48.6699, 48.9114, 49.3689, 48.6757]
    lontstart = [-93.5304, -93.5219, -93.5778, -93.0754, -92.2308]
    figfolder = fp+"/figs/20220530/"
elif (yyyyt==2022) & (mmt==5) & (ddt==21):
    lattstart = [43.0179, 42.9217, 44.1058, 44.1755]
    lontstart = [-81.2216, -81.1977, -79.1458, -78.7722]
    figfolder = fp+"/figs/20220521/"
elif (yyyyt==2025) & (mmt==7) & (ddt==24):
    lattstart = [43.4700, 44.2500, 45.0000]
    lontstart = [-81.1800, -80.5000, -79.2500]
    figfolder = fp+"/figs/20250724/"
elif (yyyyt==2026) & (mmt==4) & (ddt==15):
    lattstart = [42.2718, 42.2684, 42.2393, 42.2423]
    lontstart = [-83.7524, -83.2104, -83.0576, -82.8484]
    figfolder = fp+"/figs/20260415/"
elif (yyyyt==2026) & (mmt==6) & (ddt==30):
    lattstart = [44.4640, 44.3021, 44.3303, 45.6315, 44.6890]
    lontstart = [-76.7289, -76.5386, -76.9268, -77.8432, -76.9664]
    figfolder = fp+"/figs/20260630/"
elif (yyyyt==2026) & (mmt==7) & ((ddt==3)|(ddt==4)):
    lattstart = [42.2393, 42.4226, 42.7505, 42.2766]
    lontstart = [-83.0576, -82.1339, -81.7004, -83.7341]
    figfolder = fp+"/figs/20260703/"



[lonW,lonE,latS,latN] = [np.median(lontstart)-13, np.median(lontstart)+13, np.median(lattstart)-7, np.median(lattstart)+7]

if region=='GreatLakes':    
    [lonW,lonE,latS,latN] = [-93.0,-70.0,39.5,52.0]
    n=5
    shrinkscale = 0.6
elif region=='SasMan':
    [lonW,lonE,latS,latN] = [-114.0,-88.0,46.5,60.5]
    n=5
    shrinkscale = 0.6
elif region=='Alberta':      
    [lonW,lonE,latS,latN] = [-128.0,-102.0,46.5,60.5]
    n=5
    shrinkscale = 0.6
elif region=='Canada':
    [lonW,lonE,latS,latN] = [-132.0,-58.0,35.0,75.0]
    n=10
    shrinkscale = 0.3



timt="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))

data = xr.open_dataset(fp+f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_preslevs.nc")
datas = xr.open_dataset(fp+f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_singlevs.nc")

dat = data.sel(valid_time=timt)
dats = datas.sel(valid_time=timt)

prs = dat.pressure_level.values * units.hPa
lat = dat['latitude'].values
lon = dat['longitude'].values
X,Y = np.meshgrid(lon,lat)

data.close()
datas.close()

if hht >= 3:
    tim0="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht-3).zfill(2))
    dat0 = data.sel(valid_time=tim0)
    dats0 = datas.sel(valid_time=tim0)
else:
    tim0="%d-%s-%sT%s:00:00.000000000" %(yyyy,str(mm).zfill(2),str(dd).zfill(2),str(hht+21).zfill(2))
    data0 = xr.open_dataset(fp+f"era5_{yyyy}{mm:02.0f}{dd:02.0f}_preslevs.nc")
    datas0 = xr.open_dataset(fp+f"era5_{yyyy}{mm:02.0f}{dd:02.0f}_singlevs.nc")
    dat0 = data0.sel(valid_time=tim0)
    dats0 = datas0.sel(valid_time=tim0)
    data0.close()
    datas0.close()






for presl in [300,500,700,850]:
    
    zz = dat.sel(pressure_level=presl)['z'].values/9.81
    uu = dat.sel(pressure_level=presl)['u'].values*1.94384449
    vv = dat.sel(pressure_level=presl)['v'].values*1.94384449
    # t = dat.sel(pressure_level=presl)['t'].values-273.15
    # q = dat.sel(pressure_level=presl)['q'].values
    # omega = dat.sel(pressure_level=presl)['w'].values*units('Pa/s')
    # w = mc.vertical_velocity(omega*units('Pa/s'), np.moveaxis(np.tile(prs, (len(lat),len(lon),1)), -1, 0), t*units.degC, q*units('kg/kg'))
    
    z0 = dat0.sel(pressure_level=presl)['z'].values/9.81
    dz = zz - z0
    
    # print(f"{presl} mb")
    # print(f"           max rise: {np.max(dz):.0f} m")
    # print(f"           max fall: {np.min(dz):.0f} m")
    # print(f"           90% rise: {np.percentile(dz,90):.0f} m")
    # print(f"           10% fall: {np.percentile(dz,10):.0f} m")
    
    
    
    fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()}, layout='tight')
    ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle='-')
    ax.add_feature(cfeature.STATES, linestyle='-')
    ax.add_feature(cfeature.COASTLINE)
    
    
    colormap = pyart.graph.cmweather.cm_colorblind.balance.copy()
    colormap.set_over("maroon")
    colormap.set_under("navy")
    
    
    if presl == 300:
        levs = np.linspace(-100,100,41)
        vmin = levs[0]
        vmax = levs[-1]
        
        cf = ax.contourf(X, Y, dz, levels=levs, cmap=colormap, vmin=vmin, vmax=vmax, alpha=0.8, extend='both')
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Height change [m]', shrink=shrinkscale)
        cbar.set_label('Height change [m]', fontsize=16)  
        cbar.ax.tick_params(labelsize=12)
        
        cs = ax.contour(X, Y, zz, levels=np.arange(7500,9930,60), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
        
    
    if presl == 500:
        levs = np.linspace(-60,60,25)
        vmin = levs[0]
        vmax = levs[-1]
        
        cf = ax.contourf(X, Y, dz, levels=levs, cmap=colormap, vmin=vmin, vmax=vmax, alpha=0.8, extend='both')
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Height change [m]', shrink=shrinkscale)
        cbar.set_label('Height change [m]', fontsize=16)  
        cbar.ax.tick_params(labelsize=12)
        
        cs = ax.contour(X, Y, zz, levels=np.arange(4500,6630,60), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
        
    
    if presl == 700:
        levs = np.linspace(-40,40,17)
        vmin = levs[0]
        vmax = levs[-1]
        
        cf = ax.contourf(X, Y, dz, levels=levs, cmap=colormap, vmin=vmin, vmax=vmax, alpha=0.8, extend='both')
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Height change [m]', shrink=shrinkscale)
        cbar.set_label('Height change [m]', fontsize=16)  
        cbar.ax.tick_params(labelsize=12)
        
        cs = ax.contour(X, Y, zz, levels=np.arange(2100,3630,30), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
        
    if presl == 850:
        levs = np.linspace(-40,40,17)
        vmin = levs[0]
        vmax = levs[-1]
        
        cf = ax.contourf(X, Y, dz, levels=levs, cmap=colormap, vmin=vmin, vmax=vmax, alpha=0.8, extend='both')
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Height change [m]', shrink=shrinkscale)
        cbar.set_label('Height change [m]', fontsize=16)  
        cbar.ax.tick_params(labelsize=12)
        
        cs = ax.contour(X, Y, zz, levels=np.arange(900,1830,30), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
    
    #Wind barbs are plotted every n gridpoint
    cuv = ax.barbs(X[::n,::n], Y[::n,::n], uu[::n,::n], vv[::n,::n])
    
    # cw = ax.contour(X, Y, omega, levels=[-0.5,0.5], colors=['b','r'], linestyles='-', linewidths=1.5)
    # ax.clabel(cw, inline=True, fontsize=11, fmt='%.1f')#.1f')
    # cw2 = ax.contour(X, Y, omega, levels=[0.5,1], colors='r', linestyles='-', linewidths=1.5)
    # ax.clabel(cw2, inline=True, fontsize=11, fmt='%d')#.1f')
    
    s = ax.scatter(lontstart, lattstart, s=250, marker="v", edgecolors='k', color='yellow', linewidth=2, zorder=10)
    ax.legend(handles=[s], labels=['Tornado'], loc='upper right', fontsize=20, bbox_to_anchor=(1,1.07), frameon=False)
    
    ax.set_title(f"{presl} mb 3-h height falls, geopotential height\n {timt[:10]} {timt[11:13]}Z ", fontsize=20)
    
    if region == 'Canada':
        figname = f"heightfalls_{presl}mb_{hht:02.0f}Z_Canada.png"
    else:
        figname = f"heightfalls_{presl}mb_{hht:02.0f}Z.png"
    
    if figsave:
        plt.savefig(figfolder+figname, dpi=300)
    # plt.show()





for presl in [300,500,700,850]:
    
    zz = dat.sel(pressure_level=presl)['z'].values/9.81
    uu = dat.sel(pressure_level=presl)['u'].values*1.94384449
    vv = dat.sel(pressure_level=presl)['v'].values*1.94384449
    # t = dat.sel(pressure_level=presl)['t'].values-273.15
    # q = dat.sel(pressure_level=presl)['q'].values
    # omega = dat.sel(pressure_level=presl)['w'].values
    # w = mc.vertical_velocity(omega*units('Pa/s'), presl*100*np.ones((len(lat),len(lon)))*units.Pa, t*units.degC, q*units('kg/kg')).magnitude
    vort = dat.sel(pressure_level=presl)['vo'].values
    
    vo0 = dat0.sel(pressure_level=presl)['vo'].values
    dvort = vort - vo0
    
    # print(f"{presl} mb")
    # print(f"           max vort rise: {np.max(dvort):.7f} 1/s")
    # print(f"           min vort fall: {np.min(dvort):.7f} 1/s")
    # print(f"                max vort: {np.max(vort):.7f} 1/s")
    # print(f"                min vort: {np.min(vort):.7f} 1/s")
    # print(f"                   max w: {np.max(w):.2f} m/s")
    # print(f"                   min w: {np.min(w):.2f} m/s")
    
    
    # colormap = pyart.graph.cmweather.cm_colorblind.balance.copy()
    # colormap.set_over("maroon")
    # colormap.set_under("navy")
    colormap = cm.PuOr_r.copy()
    colormap.set_over("darkorange")
    colormap.set_under("indigo")
    
    
    fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()}, layout='tight')
    ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle='-')
    ax.add_feature(cfeature.STATES, linestyle='-')
    ax.add_feature(cfeature.COASTLINE)
    
    
    if presl == 300:
        levs = np.linspace(-6e-4,6e-4,25)
        vmin = levs[0]
        vmax = levs[-1]
        
        cf = ax.contourf(X, Y, vort, levels=levs, cmap=colormap, vmin=vmin, vmax=vmax, alpha=0.9, extend='both')
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Vorticity [1/s]', shrink=shrinkscale)
        cbar.set_label('Vorticity [1/s]', fontsize=16)
        cbar.ax.tick_params(labelsize=12)
        
        cs = ax.contour(X, Y, zz, levels=np.arange(7500,9930,60), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
        
    
    if presl == 500:
        levs = np.linspace(-4e-4,4e-4,17)
        vmin = levs[0]
        vmax = levs[-1]
        
        cf = ax.contourf(X, Y, vort, levels=levs, cmap=colormap, vmin=vmin, vmax=vmax, alpha=0.9, extend='both')
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Vorticity [1/s]', shrink=shrinkscale)
        cbar.set_label('Vorticity [1/s]', fontsize=16)
        cbar.ax.tick_params(labelsize=12)
        
        cs = ax.contour(X, Y, zz, levels=np.arange(4500,6630,60), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
        
    
    if presl == 700:
        levs = np.linspace(-3e-4,3e-4,13)
        vmin = levs[0]
        vmax = levs[-1]
        
        cf = ax.contourf(X, Y, vort, levels=levs, cmap=colormap, vmin=vmin, vmax=vmax, alpha=0.9, extend='both')
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Vorticity [1/s]', shrink=shrinkscale)
        cbar.set_label('Vorticity [1/s]', fontsize=16)
        cbar.ax.tick_params(labelsize=12)
        
        cs = ax.contour(X, Y, zz, levels=np.arange(2100,3630,30), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
        
    if presl == 850:
        levs = np.linspace(-3e-4,3e-4,13)
        vmin = levs[0]
        vmax = levs[-1]
        
        cf = ax.contourf(X, Y, vort, levels=levs, cmap=colormap, vmin=vmin, vmax=vmax, alpha=0.9, extend='both')
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Vorticity [1/s]', shrink=shrinkscale)
        cbar.set_label('Vorticity [1/s]', fontsize=16)
        cbar.ax.tick_params(labelsize=12)
        
        cs = ax.contour(X, Y, zz, levels=np.arange(900,1830,30), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
        
    cv1 = ax.contour(X, Y, dvort*1e4, levels=np.arange(-5,0), colors='b', linewidths=1.25, linestyles='-')
    ax.clabel(cv1, inline=True, fontsize=11, fmt='%d')#.1f')
    cv2 = ax.contour(X, Y, dvort*1e4, levels=np.arange(1,6), colors='r', linewidths=1.25, linestyles='-')
    ax.clabel(cv2, inline=True, fontsize=11, fmt='%d')#.1f')
    
    #Wind barbs are plotted every n gridpoint
    cuv = ax.barbs(X[::n,::n], Y[::n,::n], uu[::n,::n], vv[::n,::n])
    
    s = ax.scatter(lontstart, lattstart, s=250, marker="v", edgecolors='k', color='yellow', linewidth=2, zorder=10)
    ax.legend(handles=[s], labels=['Tornado'], loc='upper right', fontsize=20, bbox_to_anchor=(1,1.07), frameon=False)
    
    ax.set_title(f"{presl} mb vorticity, 3-H d(vorticity), geopotential height\n {timt[:10]} {timt[11:13]}Z ", fontsize=20)
    
    if region == 'Canada':
        figname = f"vorticity_{presl}mb_{hht:02.0f}Z_Canada.png"
    else:
        figname = f"vorticity_{presl}mb_{hht:02.0f}Z.png"
    
    if figsave:
        plt.savefig(figfolder+figname, dpi=300)
    # plt.show()



# data.close()
# datas.close()
# data0.close()
# datas0.close()
dat.close()
dats.close()
dat0.close()
dats0.close()

