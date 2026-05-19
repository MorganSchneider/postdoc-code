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
    [lonW,lonE,latS,latN] = [-132.0,-50.0,35.0,75.0] 
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


#%% Calculate trough tilt



for pres in [850,700,500,300]:
    dat1 = data.sel(pressure_level=pres, valid_time=timt)
    z1 = dat1["z"].values/9.81
    lon1 = dat1["longitude"].values
    lat1 = dat1["latitude"].values
    u1 = dat1["u"].values
    v1 = dat1["v"].values
    LON,LAT = np.meshgrid(lon1,lat1)
    
    xx,yy = latlon2xy(lat1, lon1, lat1[0], lon1[0])
    dvdx1 = np.zeros(shape=(len(lat1),len(lon1)))
    dudy1 = np.zeros(shape=(len(lat1),len(lon1)))
    for i in range(len(lat1)):
        dvdx1[i,:] = np.gradient(v1[i,:], xx[i,:]*1000)
    for j in range(len(lon1)):
        dudy1[:,j] = np.gradient(u1[:,j], yy[:,j]*1000)
    vort1 = dvdx1 - dudy1
    # X,Y = np.meshgrid(xx,yy)
    
    vortmin = np.where(vort1 == np.min(vort1))
    ixmin = vortmin[1][0]
    iymin = vortmin[0][0]
    vortmax = np.where(vort1 == np.max(vort1))
    ixmax = vortmax[1][0]
    iymax = vortmax[0][0]
    
    xmin = xx[iymin,ixmin]
    ymin = yy[iymin,ixmin]
    xmax = xx[iymax,ixmax]
    ymax = yy[iymax,ixmax]
    # xmin,ymin = latlon2xy(lat1[iymin], lon1[ixmin], lat1[0], lon1[0])
    # xmax,ymax = latlon2xy(lat1[iymax], lon1[ixmax], lat1[0], lon1[0])
    
    if pres == 850:
        trough_tilt_850 = (ymax-ymin) / (xmax-xmin)
    elif pres == 700:
        trough_tilt_700 = (ymax-ymin) / (xmax-xmin)
    elif pres == 500:
        trough_tilt_500 = (ymax-ymin) / (xmax-xmin)
    elif pres == 300:
        trough_tilt_300 = (ymax-ymin) / (xmax-xmin)
    
    # izmin = np.where(z1 == np.min(z1))
    # iy0 = izmin[0][0]
    # ix0 = izmin[1][0]
    # x_trough = np.zeros(shape=(iy0,))
    # y_trough = np.zeros(shape=(iy0,))
    # for ii in range(iy0):
    #     ix = np.argmin(abs(z1[ii,:] - np.min(z1[ii,:])))
    #     # x,y = latlon2xy(lat1[ii], lon1[ix], lat1[0], lon1[0])
    #     x_trough[ii] = xx[ix,ii]
    #     y_trough[ii] = yy[ix,ii]
    
    # slopes = np.zeros(shape=(len(x_trough),))
    # for i in range(len(x_trough)-1):
    #     slopes[i] = (y_trough[i+1]-y_trough[i]) / (x_trough[i+1]-x_trough[i])
    
    # if pres == 850:
    #     trough_tilt_850 = np.mean(slopes)
    # elif pres == 700:
    #     trough_tilt_700 = np.mean(slopes)
    # elif pres == 500:
    #     trough_tilt_500 = np.mean(slopes)
    # elif pres == 300:
    #     trough_tilt_300 = np.mean(slopes)







