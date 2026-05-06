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

case='tornado'
#case='hail'

region='Alberta'
#region='SasMan'
#region='GreatLakes'
#region='Canada'

yyyyt = 2023
mmt   = 7
ddt   = 1
hht   = 20
lattstart=51.6107  # from NTP/NHP website
lontstart=-114.198 # from NTP/NHP website

timt="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))

lattstart=51.6107  #Latitude of tornado/beginning of hail swath (from NTP/NHP database)
lontstart=-114.198 #Longitude of tornado/beginning of hail swath (from NTP/NHP database)

latt=52.25       # latitude of sounding
lont=-113.25     # longitude of sounding #if you do not need this, uncomment the next two lines
# latt=lattstart 
# lont=lontstart

figfolder='./figures/'
data_preslevs=f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_preslevs.nc"
data_singlevs=f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_singlevs.nc"


#%% Define plot domain - can add other regions here

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
    [lonW,lonE,latS,latN] = [-132.0,-50.0,40.0,65.0] 
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

PlotCAPEMap(data,datas,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,figfolder)

PlotPressureMaps(data,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,figfolder)

PlotMoistureMap2(data,datas,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,figfolder)

data00 = data.sel(valid_time=timt)
z00=data00['z'].values/9.81
u00=data00['u'].values
v00=data00['v'].values
lat=data00['latitude'].values
lon=data00['longitude'].values

latt=51.5  #Latitude of tornado/beginning of hail swath (from NTP/NHP database)
lont=-114.0

lati=np.where(lat==latt)
loni=np.where(lon==lont)

print(np.shape(z00))
print(lati,loni)

for i in range(0,np.shape(z00)[0]):
    print('%d\t%.1f\t%.1f' %(z00[i,lati,loni],u00[i,lati,loni],v00[i,lati,loni]))










