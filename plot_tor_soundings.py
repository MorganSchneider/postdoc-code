# -*- coding: utf-8 -*-
"""
Created on Wed May 13 08:32:35 2026

@author: mschne28
"""

# from era5utils import *

from CM1utils import *

import xarray as xr
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, BoundaryNorm
from matplotlib import  cm
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors

import metpy.calc as mc
from metpy.cbook import get_test_data
from metpy.plots import Hodograph, SkewT
from metpy.units import units
# from metpy.calc import dewpoint_from_relative_humidity,specific_humidity_from_dewpoint
# from metpy.calc import potential_temperature,temperature_from_potential_temperature
# from metpy.calc import wind_speed, wind_direction,bunkers_storm_motion

from era5utils import extract_data, plot_with_hodograph, plot_skewT



#%% Load data and plot one sounding


fp = "C:/Users/mschne28/Documents/era5/tor_outbreaks/"

yyyy = 2021
mm = 8
dd = 11

tor_num = 1


dbfile = open(fp+"tornado_locs.pkl", 'rb')
locs_all = pickle.load(dbfile)
locs = locs_all[f"{yyyy}{mm:02.0f}{dd:02.0f}"]
dbfile.close()

tloc = locs[f"loc{tor_num}"]

name = tloc["name"]
yyyyt = tloc["date_ymd"][0]
mmt = tloc["date_ymd"][1]
ddt = tloc["date_ymd"][2]
timestr = tloc["time_utc"]
lattstart = tloc["lat"]
lontstart = tloc["lon"]

hht = float(timestr) // 100
# latt = np.round(lattstart/0.25) * 0.25
# lont = np.round(lontstart/0.25) * 0.25


fn_preslev = f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_preslevs.nc"
fn_singlev = f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_singlevs.nc"

# timt="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))
timt = f"{yyyyt}-{mmt:02.0f}-{ddt:02.0f}T{hht:02.0f}:00:00.000000000"


dsp = xr.open_dataset(fp+fn_preslev)
dss = xr.open_dataset(fp+fn_singlev)

lats = dsp['latitude'][:].values
lons = dsp['longitude'][:].values

lati = np.argmin(np.abs(lats-lattstart))
loni = np.argmin(np.abs(lons-lontstart))

latt = lats[lati]
lont = lons[loni]

# p,z,T,q,theta,Td,u,v,u10,v10,speed,direc,cape,cin,sfcp,orog,q2m,theta2m,td2m,t2m,leftm,meanm,rightm,parcel_prof,lcl_pressure,lcl_temperature = extract_data(latt,lont,timt,dsp,dss)

data = extract_data(latt, lont, timt, dsp, dss)
# z = data['z']
# orog = data['orog']
# T = data['T']
# Td = data['Td']
# p = data['p']
# u = data['u']
# v = data['v']
# u10 = data['u10']
# v10 = data['v10']
# leftmover = data['leftm']
# meanwind = data['meanm']
# rightmover = data['rightm']
# lcl_pressure = data['lcl_pressure']
# lcl_temperature = data['lcl_temperature']
# parcel_prof = data['parcel_prof']
figfolder = fp+"figs/"

dsp.close()
dss.close()




figsave = False


fig = plot_skewT(timt, latt, lont, data, tloc, figfolder=None, figsave=False)


# plot_with_hodograph(timt, lont, latt, z, orog, T, Td, p, u, v, u10, v10, leftm, meanm, rightm, lcl_pressure, lcl_temperature, parcel_prof, pngfolder=figfolder)

# fig = plt.figure(figsize=(9,9))
# skew = SkewT(fig)

# skew.ax.set_ylim(1000, 100)
# skew.ax.set_xlim(-50, 40)

# # Plot special lines
# skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black') # Plot LCL as black dot
# skew.ax.axvline(0, color='c', linestyle='--', linewidth=2) # Plot a zero degree isotherm
# skew.plot_dry_adiabats()
# skew.plot_moist_adiabats()
# skew.plot_mixing_lines()

# # Plot the data using normal plotting functions, in this case using
# # log scaling in Y, as dictated by the typical meteorological plot
# skew.plot(p[z>orog], T[z>orog], 'r', linewidth=2)
# skew.plot(p[z>orog], Td[z>orog], 'g', linewidth=2)
# skew.plot(p[z>orog], parcel_prof, 'k', linewidth=2) # Plot the parcel profile as a black line
# skew.shade_cin(p[z>orog], T[z>orog], parcel_prof, Td[z>orog]) # Shade areas of CAPE and CIN
# skew.shade_cape(p[z>orog], T[z>orog], parcel_prof)
# skew.plot_barbs(p[z>orog], u[z>orog], v[z>orog], xloc=1.1, plot_units=units('kts'))


# skew.ax.set_xlabel("Temperature (C)", fontsize=12)
# skew.ax.set_ylabel("Pressure (hPa)", fontsize=12)
# # skew.ax.set_title(f"{name} , {timt[:13]} UTC \n lat,lon = {latt:.2f}, {lont:.2f} ", fontsize=12)
# skew.ax.set_title(f"{name} , {yyyyt}-{mmt:02.0f}-{ddt:02.0f} T{timestr[:2]}:{timestr[2:]} UTC ({lattstart:.2f}, {lontstart:.2f}) \n ERA5 sounding: {timt[:13]} UTC ({latt:.2f}, {lont:.2f})", fontsize=12)

# # Create a hodograph
# ax_hod = inset_axes(skew.ax, '40%', '40%', loc=1)
# h = Hodograph(ax_hod, component_range=40.)
# h.add_grid(increment=10)

# cmap = plt.get_cmap('tab20') # Get the tab20 colormap from matplotlib
# colors = [cmap(i) for i in range(7)] # Extract the first 7 colors from it
# my_cmap = ListedColormap(colors, name="my_cmap")
# my_cmap_r = my_cmap.reversed()

# unew = np.zeros((len(u[z>orog])+1))
# vnew = np.zeros((len(v[z>orog])+1))
# znew = np.zeros((len(z[z>orog])+1))
# for k in range(0,len(unew)):
#     if k == 0:
#         unew[k] = u10.magnitude
#         vnew[k] = v10.magnitude
#         znew[k] = 10
#     else:
#         unew[k] = u[z>orog][k-1].magnitude
#         vnew[k] = v[z>orog][k-1].magnitude
#         znew[k] = z[z>orog][k-1] - orog

# hplot = h.plot_colormapped(unew*units("m/s"), vnew*units("m/s"), znew/1000, intervals=[0,1,3,5,7.5,10,12.5], colors=colors)

# cax = ax_hod.inset_axes([0.06, 0.1, 0.88, 0.04])    
# fig.colorbar(hplot, cax=cax, orientation='horizontal',ticks=[0,1,3,5,7.5,10,12.5])

# ax_hod.scatter(leftmover.magnitude[0], leftmover.magnitude[1], s=10, color='k')
# ax_hod.scatter(rightmover.magnitude[0], rightmover.magnitude[1], s=10, color='k')
# ax_hod.scatter(meanwind.magnitude[0], meanwind.magnitude[1], s=10, color='k')
# # ax_hod.scatter(u10.magnitude, v10.magnitude, s=25, color='k')

# ax_hod.text(leftmover.magnitude[0], leftmover.magnitude[1], 'LM ', ha='right', va='center_baseline', fontsize=8, fontweight='bold')
# ax_hod.text(rightmover.magnitude[0], rightmover.magnitude[1], 'RM ', ha='right', va='center_baseline', fontsize=8, fontweight='bold')
# ax_hod.text(meanwind.magnitude[0], meanwind.magnitude[1], 'MW ', ha='right', va='center_baseline', fontsize=8, fontweight='bold')

# ax_hod.xaxis.set_major_locator(MultipleLocator(20))
# ax_hod.xaxis.set_minor_locator(MultipleLocator(10))
# ax_hod.yaxis.set_major_locator(MultipleLocator(20))
# ax_hod.yaxis.set_minor_locator(MultipleLocator(10))

# if figsave:
#     plt.savefig(fp+f"figs/input_sounding_{timt[:13]}_{latt:.2f}_{lont:.2f}.png", facecolor='white', transparent=False)
#     # plt.close()

# plt.show()

#%% Load data and plot all soundings per event

fp = "C:/Users/mschne28/Documents/era5/tor_outbreaks/"

yyyy = 2021
mm = 8
dd = 11

ntors = 6


dbfile = open(fp+"tornado_locs.pkl", 'rb')
locs_all = pickle.load(dbfile)
locs = locs_all[f"{yyyy}{mm:02.0f}{dd:02.0f}"]
dbfile.close()

for i in range(ntors):
    tloc = locs[f"loc{i+1}"]
    
    name = tloc["name"]
    yyyyt = tloc["date_ymd"][0]
    mmt = tloc["date_ymd"][1]
    ddt = tloc["date_ymd"][2]
    timestr = tloc["time_utc"]
    lattstart = tloc["lat"]
    lontstart = tloc["lon"]
    
    hht = float(timestr) // 100
    # latt = np.round(lattstart/0.25) * 0.25
    # lont = np.round(lontstart/0.25) * 0.25


    fn_preslev = f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_preslevs.nc"
    fn_singlev = f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_singlevs.nc"

    # timt="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))
    timt = f"{yyyyt}-{mmt:02.0f}-{ddt:02.0f}T{hht:02.0f}:00:00.000000000"


    dsp = xr.open_dataset(fp+fn_preslev)
    dss = xr.open_dataset(fp+fn_singlev)

    lats = dsp['latitude'][:].values
    lons = dsp['longitude'][:].values

    lati = np.argmin(np.abs(lats-lattstart))
    loni = np.argmin(np.abs(lons-lontstart))

    latt = lats[lati]
    lont = lons[loni]

    # p,z,T,q,theta,Td,u,v,u10,v10,speed,direc,cape,cin,sfcp,orog,q2m,theta2m,td2m,t2m,leftm,meanm,rightm,parcel_prof,lcl_pressure,lcl_temperature = extract_data(latt,lont,timt,dsp,dss)

    data = extract_data(latt, lont, timt, dsp, dss)
    # z = data['z']
    # orog = data['orog']
    # T = data['T']
    # Td = data['Td']
    # p = data['p']
    # u = data['u']
    # v = data['v']
    # u10 = data['u10']
    # v10 = data['v10']
    # leftmover = data['leftm']
    # meanwind = data['meanm']
    # rightmover = data['rightm']
    # lcl_pressure = data['lcl_pressure']
    # lcl_temperature = data['lcl_temperature']
    # parcel_prof = data['parcel_prof']
    figfolder = None

    dsp.close()
    dss.close()
    
    figsave = False
    
    
    fig = plot_skewT(timt, latt, lont, data, tloc, figfolder=None, figsave=False)
    
    
    
    































