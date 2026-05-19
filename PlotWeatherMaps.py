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
# region='GreatLakes'
region='Canada'

yyyyt = 2021
mmt   = 8
ddt   = 11
hht   = 20

timt="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))


lattstart = 46.4797  #Latitude of tornado/beginning of hail swath (from NTP/NHP database)
lontstart = -83.6403 #Longitude of tornado/beginning of hail swath (from NTP/NHP database)


latt = 46.50       # latitude of sounding
lont = -83.75     # longitude of sounding #if you do not need this, uncomment the next two lines
# latt=lattstart 
# lont=lontstart

figfolder=fp+'figs/'
data_preslevs=fp+f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_preslevs.nc"
data_singlevs=fp+f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_singlevs.nc"


#% Define plot domain - can add other regions here

#Plotdomain:
#lonW=-125
#lonE=-90
#latS=45
#latN=65

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



#%% Read data files

# Pressure level data
data = xr.open_dataset(data_preslevs)
data.head(5)

# Single level data
datas = xr.open_dataset(data_singlevs)
datas.head(5)


#%% Make plots

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


#%% Calculate trough tilt?

p850 = {'latmin':0, 'lonmin':0, 'latmax':0, 'lonmax':0, 'xmin':0, 'ymin':0, 'xmax':0, 'ymax':0, 'tilt':0}
p700 = {'latmin':0, 'lonmin':0, 'latmax':0, 'lonmax':0, 'xmin':0, 'ymin':0, 'xmax':0, 'ymax':0, 'tilt':0}
p500 = {'latmin':0, 'lonmin':0, 'latmax':0, 'lonmax':0, 'xmin':0, 'ymin':0, 'xmax':0, 'ymax':0, 'tilt':0}
p300 = {'latmin':0, 'lonmin':0, 'latmax':0, 'lonmax':0, 'xmin':0, 'ymin':0, 'xmax':0, 'ymax':0, 'tilt':0}
trough = {'850':{}, '700':{}, '500':{}, '300':{}}



for pres in [850,700,500,300]:
    dat1 = data.sel(pressure_level=pres, valid_time=timt)
    z1 = dat1["z"].values/9.81
    lon1 = dat1["longitude"].values
    lat1 = dat1["latitude"].values
    u1 = dat1["u"].values
    v1 = dat1["v"].values
    vort1 = dat1["vo"].values
    pv1 = dat1["pv"].values
    LON,LAT = np.meshgrid(lon1,lat1)
    
    # xx,yy = latlon2xy(lat1, lon1, lat1[0], lon1[0])
    # X,Y = np.meshgrid(xx,yy)
    
    imin = np.where((z1 == np.min(z1[(lat1<=65),:])) & (LAT<=65))
    latmin = lat1[imin[0][0]]
    lonmin = lon1[imin[1][0]]
    
    imax = np.where((z1 == np.max(z1[:,(lon1>lonmin)])) & (LON>lonmin))
    latmax = lat1[imax[0][0]]
    lonmax = lon1[imax[1][0]]
    
    
    # xmin = xx[iymin,ixmin]
    # ymin = yy[iymin,ixmin]
    # xmax = xx[iymax,ixmax]
    # ymax = yy[iymax,ixmax]
    xmin,ymin = latlon2xy(latmin, lonmin, lat1[-1], lon1[0])
    xmax,ymax = latlon2xy(latmax, lonmax, lat1[-1], lon1[0])
    
    trough_tilt = (ymax-ymin) / (xmax-xmin)
    
    trough[f"{pres}"].update({'lonmin':lonmin, 'latmin':latmin, 'latmax':latmax, 'lonmax':lonmax,
                              'xmin':xmin, 'ymin':ymin, 'xmax':xmax, 'ymax':ymax,
                              'tilt':trough_tilt})
    

print(f"Trough tilt (850 mb) = {trough['850']['tilt']}")
print(f"Trough tilt (700 mb) = {trough['700']['tilt']}")
print(f"Trough tilt (500 mb) = {trough['500']['tilt']}")
print(f"Trough tilt (300 mb) = {trough['300']['tilt']}")


