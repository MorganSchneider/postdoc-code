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

from era5utils import *



#%% Load data and plot one sounding


fp = "C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/tor_outbreaks/"

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
figfolder = fp+f"figs/{yyyy}{mm:02.0f}{dd:02.0f}/"

# dsp.close()
# dss.close()




# figsave = False


# fig = plot_skewT(timt, latt, lont, data, tloc, figfolder=None, figsave=False)


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

fp = "C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/tor_outbreaks/"
# figfolder = fp+f"figs/"

yyyy = 2021; mm = 8; dd = 11; ntors = 6; casetype = "Tornado outbreak"
# yyyy = 2025; mm = 6; dd = 23; ntors = 6; casetype = "Tornado outbreak"
# yyyy = 2022; mm = 5; dd = 30; ntors = 5; casetype = "Tornado sub-outbreak"
# yyyy = 2022; mm = 5; dd = 21; ntors = 4; casetype = "Tornado sub-outbreak"
# yyyy = 2026; mm = 6; dd = 30; ntors = 5; casetype = "Tornado sub-outbreak"
# yyyy = 2025; mm = 7; dd = 24; ntors = 3; casetype = "Null event"
# yyyy = 2026; mm = 7; dd = 3; ntors = 4; casetype = "Null event"
# yyyy = 2026; mm = 4; dd = 15; ntors = 4; casetype = "Null event"


# yyyyt = 2021; mmt = 8; ddt = 11; tor hht = 20,21 (17-22)
# yyyyt = 2025; mmt = 6; ddt = 23-24; tor hht = 20,22,1 (17-23, 0-3)
# yyyyt = 2022; mmt = 5; ddt = 30-31; tor hht = 0,1,2 (21-23, 0-3)
# yyyyt = 2022; mmt = 5; ddt = 21; tor hht = 15,17 (12-18)
# yyyyt = 2026; mmt = 6; ddt = 30; tor hht = 16,17 (13-17)
# yyyyt = 2025; mmt = 7; ddt = 24-25; bow/db hht = 21-23 (19-23, 0-1)
# yyyyt = 2026; mmt = 7; ddt = 3-4; hht = 22-00 (19-23, 0-2)
# yyyyt = 2026; mmt = 4; ddt = 15; tor/null hht = 5,6 (2-8)



figfolder = fp+f"figs/{yyyy}{mm:02.0f}{dd:02.0f}/"

leadtime = 1 #hours before hit


# ntors = 1



figsave = False



dbfile = open(fp+"tornado_locs.pkl", 'rb')
locs_all = pickle.load(dbfile)
locs = locs_all[f"{yyyy}{mm:02.0f}{dd:02.0f}"]
dbfile.close()

z = np.zeros((37,ntors))
orog = np.zeros((ntors,))
T = np.zeros((37,ntors))
Td = np.zeros((37,ntors))
u = np.zeros((37,ntors))
v = np.zeros((37,ntors))
u10 = np.zeros((ntors,))
v10 = np.zeros((ntors,))
leftmover = np.zeros((2,ntors))
meanwind = np.zeros((2,ntors))
rightmover = np.zeros((2,ntors))
cape = np.zeros((ntors,))
cin = np.zeros((ntors,))
lcl_pressure = np.zeros((ntors,))
lcl_temperature = np.zeros((ntors,))
parcel_prof = np.zeros((37,ntors))

q = np.zeros((37,ntors))
theta = np.zeros((37,ntors))
speed = np.zeros((37,ntors))
direc = np.zeros((37,ntors))
sfcp = np.zeros((ntors,))
q2m = np.zeros((ntors,))
t2m = np.zeros((ntors,))
td2m = np.zeros((ntors,))
theta2m = np.zeros((ntors,))

shear01 = np.zeros((ntors,))
shear03 = np.zeros((ntors,))
shear06 = np.zeros((ntors,))
srh01 = np.zeros((ntors,))
srh03 = np.zeros((ntors,))
lapse_rate = np.zeros((ntors,))



for i in range(ntors):
    tloc = locs[f"loc{i+1}"]
    
    # tloc = locs[f"loc{i+3}"]
    
    name = tloc["name"]
    yyyyt = tloc["date_ymd"][0]
    mmt = tloc["date_ymd"][1]
    ddt = tloc["date_ymd"][2]
    timestr = tloc["time_utc"]
    lattstart = tloc["lat"]
    lontstart = tloc["lon"]
    
    hht = float(timestr) // 100 - leadtime
    # latt = np.round(lattstart/0.25) * 0.25
    # lont = np.round(lontstart/0.25) * 0.25


    fn_preslev = f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_preslevs.nc"
    fn_singlev = f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_singlevs.nc"

    # timt="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))
    timt = f"{yyyyt}-{mmt:02.0f}-{ddt:02.0f}T{hht:02.0f}:00:00.000000000"


    datap = xr.open_dataset(fp+fn_preslev)
    datas = xr.open_dataset(fp+fn_singlev)

    latitude = datap['latitude'][:].values
    longitude = datap['longitude'][:].values

    lati = np.argmin(np.abs(latitude-lattstart))
    loni = np.argmin(np.abs(longitude-lontstart))
    
    # latt = latitude[lati]
    # lont = longitude[loni]
    latt = latitude[lati-4:lati+4]
    lont = longitude[loni-4:loni+4]
    
    # p,z,T,q,theta,Td,u,v,u10,v10,speed,direc,cape,cin,sfcp,orog,q2m,theta2m,td2m,t2m,leftm,meanm,rightm,parcel_prof,lcl_pressure,lcl_temperature = extract_data(latt,lont,timt,dsp,dss)
    data = extract_data(latt, lont, timt, datap, datas)
    
    z[:,i] = data['z']
    orog[i] = data['orog']
    T[:,i] = data['T']
    Td[:,i] = data['Td']
    p = data['p']
    u[:,i] = data['u']
    v[:,i] = data['v']
    u10[i] = data['u10'].magnitude
    v10[i] = data['v10'].magnitude
    leftmover[:,i] = data['leftm']
    meanwind[:,i] = data['meanm']
    rightmover[:,i] = data['rightm']
    cape[i] = data['cape']
    cin[i] = data['cin']
    lcl_pressure[i] = data['lcl_pressure'].magnitude
    lcl_temperature[i] = data['lcl_temperature'].magnitude
    parcel_prof[:,i] = data['parcel_prof'].magnitude
    q[:,i] = data['q']
    theta[:,i] = data['theta']
    speed[:,i] = data['speed']
    direc[:,i] = data['direc']
    sfcp[i] = data['sfcp']
    q2m[i] = data['q2m'].magnitude
    t2m[i] = data['t2m'].magnitude
    td2m[i] = data['td2m'].magnitude
    theta2m[i] = data['theta2m'].magnitude
    shear06[i] = data['shear06'].magnitude
    shear03[i] = data['shear03'].magnitude
    shear01[i] = data['shear01'].magnitude
    srh03[i] = data['srh03'].magnitude
    srh01[i] = data['srh01'].magnitude
    lapse_rate[i] = data['lapse_rate']
    

    datap.close()
    datas.close()
    
    
    
    # fig = plot_skewT(timt, latitude[lati], longitude[loni], data, tloc, figfolder=figfolder, figsave=figsave)



composite_data = {'p':p, 'z':np.mean(z,axis=1), 'T':np.mean(T,axis=1)*units.K, 'q':np.mean(q,axis=1), 'theta':np.mean(theta,axis=1)*units.K, 'Td':np.mean(Td,axis=1)*units.degC,
            'u':np.mean(u,axis=1)*units("m/s"), 'v':np.mean(v,axis=1)*units("m/s"), 'u10':np.mean(u10)*units("m/s"), 'v10':np.mean(v10)*units("m/s"),
            'speed':np.mean(speed,axis=1)*units("m/s"), 'direc':np.mean(direc,axis=1)*units.degree,
            'cape':np.mean(cape), 'cin':np.mean(cin), 'sfcp':np.mean(sfcp), 'orog':np.mean(orog),
            'q2m':np.mean(q2m), 'theta2m':np.mean(theta2m)*units.K, 'td2m':np.mean(td2m)*units.K, 't2m':np.mean(t2m)*units.K,
            'leftm':np.mean(leftmover,axis=1)*units("m/s"), 'meanm':np.mean(meanwind,axis=1)*units("m/s"), 'rightm':np.mean(rightmover,axis=1)*units("m/s"),
            'parcel_prof':np.nanmean(parcel_prof,axis=1)*units.degC, 'lcl_pressure':np.mean(lcl_pressure)*units.hPa, 'lcl_temperature':np.mean(lcl_temperature)*units.K,
            'shear01':np.nanmean(shear01)*units('m/s'), 'shear03':np.nanmean(shear03)*units('m/s'), 'shear06':np.nanmean(shear06)*units('m/s'),
            'srh01':np.nanmean(srh01)*units('m**2 / s**2'), 'srh03':np.nanmean(srh03)*units('m**2 / s**2'), 'lapse_rate':np.nanmean(lapse_rate)}




figsave = False

titlestr = f"{yyyy}-{mm:02.0f}-{dd:02.0f} {casetype} \n Composite ERA5 sounding, lead time T-{leadtime}H"
figname = f"composite_sounding_{yyyy}-{mm:02.0f}-{dd:02.0f}T-{leadtime}H"
fig = plot_skewT(timt, latitude[lati], longitude[loni], composite_data, tloc, titlestr=titlestr, figname=figname, figfolder=figfolder, figsave=figsave)

#%% skew t
z = composite_data['z']
orog = composite_data['orog']
p = composite_data['p']
T = composite_data['T']
Td = composite_data['Td']
parcel_prof = composite_data['parcel_prof']
u = composite_data['u']
v = composite_data['v']
u10 = composite_data['u10']
v10 = composite_data['v10']
lcl_pressure = composite_data['lcl_pressure']
lcl_temperature = composite_data['lcl_temperature']
cape = composite_data['cape']
cin = composite_data['cin']
shear01 = composite_data['shear01']
shear03 = composite_data['shear03']
shear06 = composite_data['shear06']
srh01 = composite_data['srh01']
srh03 = composite_data['srh03']
leftmover = composite_data['leftm']
rightmover = composite_data['rightm']
meanwind = composite_data['meanm']
lm_speed = mc.wind_speed(leftmover[0], leftmover[1])
rm_speed = mc.wind_speed(rightmover[0], rightmover[1])
mw_speed = mc.wind_speed(meanwind[0], meanwind[1])



fig = plt.figure(figsize=(9,9))
skew = SkewT(fig)

skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 40)

# Plot special lines
skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black') # Plot LCL as black dot
# skew.ax.axvline(0, color='c', linestyle='--', linewidth=1.25) # Plot a zero degree isotherm
skew.plot_dry_adiabats(linewidth=1.25)
skew.plot_moist_adiabats(linewidth=1.25)
skew.plot_mixing_lines(linewidth=1.25)

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(p[z>orog], T[z>orog], 'r', linewidth=2)
skew.plot(p[z>orog], Td[z>orog], 'g', linewidth=2)
skew.plot(p[z>orog], parcel_prof[z>orog], 'k', linewidth=2) # Plot the parcel profile as a black line
skew.shade_cin(p[z>orog], T[z>orog], parcel_prof[z>orog], Td[z>orog]) # Shade areas of CAPE and CIN
skew.shade_cape(p[z>orog], T[z>orog], parcel_prof[z>orog])
skew.plot_barbs(p[z>orog], u[z>orog], v[z>orog], xloc=1.1, plot_units=units('kts'))


skew.ax.set_xlabel("Temperature (C)", fontsize=12)
skew.ax.set_ylabel("Pressure (hPa)", fontsize=12)
# skew.ax.set_title(f"{name} , {timt[:13]} UTC \n lat,lon = {latt:.2f}, {lont:.2f} ", fontsize=12)
skew.ax.set_title(titlestr, fontsize=12)

# Create a hodograph
ax_hod = inset_axes(skew.ax, '40%', '40%', loc=1)
h = Hodograph(ax_hod, component_range=40.)
h.add_grid(increment=10)

cmap = plt.get_cmap('tab20') # Get the tab20 colormap from matplotlib
colors = [cmap(i) for i in range(7)] # Extract the first 7 colors from it
my_cmap = ListedColormap(colors, name="my_cmap")
my_cmap_r = my_cmap.reversed()

unew = np.zeros((len(u)+1))
vnew = np.zeros((len(v)+1))
znew = np.zeros((len(z)+1))
for k in range(0,len(unew)):
    if k == 0:
        unew[k] = u10.magnitude
        vnew[k] = v10.magnitude
        znew[k] = 10
    else:
        unew[k] = u[k-1].magnitude
        vnew[k] = v[k-1].magnitude
        znew[k] = z[k-1] - orog

hplot = h.plot_colormapped(unew[znew>orog]*units("m/s"), vnew[znew>orog]*units("m/s"), znew[znew>orog]/1000, intervals=[0,1,3,5,7.5,10,12.5], colors=colors)

cax = ax_hod.inset_axes([0.06, 0.1, 0.88, 0.04])    
fig.colorbar(hplot, cax=cax, orientation='horizontal',ticks=[0,1,3,5,7.5,10,12.5])

ax_hod.scatter(leftmover.magnitude[0], leftmover.magnitude[1], s=10, color='k')
ax_hod.scatter(rightmover.magnitude[0], rightmover.magnitude[1], s=10, color='k')
ax_hod.scatter(meanwind.magnitude[0], meanwind.magnitude[1], s=10, color='k')
# ax_hod.scatter(u10.magnitude, v10.magnitude, s=25, color='k')

ax_hod.text(leftmover.magnitude[0], leftmover.magnitude[1], 'LM ', ha='right', va='center_baseline', fontsize=8, fontweight='bold')
ax_hod.text(rightmover.magnitude[0], rightmover.magnitude[1], 'RM ', ha='right', va='center_baseline', fontsize=8, fontweight='bold')
ax_hod.text(meanwind.magnitude[0], meanwind.magnitude[1], 'MW ', ha='right', va='center_baseline', fontsize=8, fontweight='bold')

ax_hod.xaxis.set_major_locator(MultipleLocator(20))
ax_hod.xaxis.set_minor_locator(MultipleLocator(10))
ax_hod.yaxis.set_major_locator(MultipleLocator(20))
ax_hod.yaxis.set_minor_locator(MultipleLocator(10))


str1 = '\n'.join((
    "CAPE:               %.0f J/kg" % (cape, ),
    "CIN:                 -%.0f J/kg" % (cin, ),
    "0-6km MW:       %.1f m/s" % (mw_speed.magnitude,),
    "Bunkers RM:     %.1f m/s" % (rm_speed.magnitude,),
    "Shear 0-1 km:   %.1f m/s" % (shear01.magnitude,),
    "Shear 0-3 km:   %.1f m/s" % (shear03.magnitude,),
    "Shear 0-6 km:   %.1f m/s" % (shear06.magnitude,),
    "SRH 0-1 km:     %.0f m2/s2" % (srh01.magnitude,),
    "SRH 0-3 km:     %.0f m2/s2" % (srh03.magnitude,),
    "LCL pressure:   %.0f hPa" % (lcl_pressure.magnitude,)
    ))
props = dict(boxstyle='square', facecolor='w', edgecolor='k')
skew.ax.text(-39, 950, str1, fontsize=11, bbox=props)


if figsave:
    plt.savefig(figfolder+figname+'.png', facecolor='white', transparent=False, dpi=300)

plt.show()


# print(f"---{mm:02.0f}-{dd:02.0f}-{yyyy}---")
# print(f"   Lead time: {leadtime} h")
# print(f"CAPE = {cape:.1f} J/kg")
# print(f"CIN = -{cin:.1f} J/kg")
# print(f"Shear 0-1km = {shear01.magnitude:.1f} m/s")
# print(f"Shear 0-3km = {shear03.magnitude:.1f} m/s")
# print(f"Shear 0-6km = {shear06.magnitude:.1f} m/s")
# print(f"SRH 0-3km = {srh03.magnitude:.1f} m2/s2")
# print(f"SRH 0-1km = {srh01.magnitude:.1f} m2/s2")
# print(f"LCL pres = {lcl_pressure.magnitude:.1f} hPa")







#%% Get sounding parameters - first tornado only

fp = "C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/tor_outbreaks/"


dbfile = open(fp+"tornado_locs.pkl", 'rb')
locs = pickle.load(dbfile)
# locs = locs_all[f"{yyyy}{mm:02.0f}{dd:02.0f}"]
dbfile.close()

events = ["20210811", "20250623", "20220530", "20220521", "20260630", "20250724", "20260703"]

t2m = np.zeros((len(events),))
td2m = np.zeros((len(events),))
cape = np.zeros((len(events),))
cin = np.zeros((len(events),))
shear01 = np.zeros((len(events),))
shear03 = np.zeros((len(events),))
shear06 = np.zeros((len(events),))
srh01 = np.zeros((len(events),))
srh03 = np.zeros((len(events),))
lclp = np.zeros((len(events),))
lr = np.zeros((len(events),))
sfcp = np.zeros((len(events),))
lclz = np.zeros((len(events),))
dcape = np.zeros((len(events),))
downT = np.zeros((len(events),)) #surface temperature of downdraft parcel - cold pool temperature


leadtime = 1
from era5utils import *

for i in range(len(events)):
    tor1 = locs[events[i]]['loc1']
    
    name = tor1['name']
    yyyyt = tor1["date_ymd"][0]
    mmt = tor1["date_ymd"][1]
    ddt = tor1["date_ymd"][2]
    timestr = tor1["time_utc"]
    lattstart = tor1["lat"]
    lontstart = tor1["lon"]
    hht = float(timestr) // 100 - leadtime
    
    fn_preslev = f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_preslevs.nc"
    fn_singlev = f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_singlevs.nc"

    # timt="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))
    timt = f"{yyyyt}-{mmt:02.0f}-{ddt:02.0f}T{hht:02.0f}:00:00.000000000"
    
    datap = xr.open_dataset(fp+fn_preslev)
    datas = xr.open_dataset(fp+fn_singlev)
    
    latitude = datap['latitude'][:].values
    longitude = datap['longitude'][:].values
    
    lati = np.argmin(np.abs(latitude-lattstart))
    loni = np.argmin(np.abs(longitude-lontstart))
    
    latt = latitude[lati-3:lati+4]
    lont = longitude[loni-3:loni+4]
    
    # p,z,T,q,theta,Td,u,v,u10,v10,speed,direc,cape,cin,sfcp,orog,q2m,theta2m,td2m,t2m,leftm,meanm,rightm,parcel_prof,lcl_pressure,lcl_temperature = extract_data(latt,lont,timt,dsp,dss)
    data = extract_data(latt, lont, timt, datap, datas)
    
    p = data['p'].magnitude
    t2m[i] = data['t2m'].magnitude
    td2m[i] = data['td2m'].magnitude
    cape[i] = data['cape']
    cin[i] = data['cin']
    lclp[i] = data['lcl_pressure'].magnitude
    shear06[i] = data['shear06'].magnitude
    shear03[i] = data['shear03'].magnitude
    shear01[i] = data['shear01'].magnitude
    srh03[i] = data['srh03'].magnitude
    srh01[i] = data['srh01'].magnitude
    lr[i] = data['lapse_rate']
    sfcp[i] = data['sfcp']
    lclz[i] = data['lcl_height'].magnitude
    dcape[i] = data['dcape']
    downT[i] = data['downT'].magnitude
        
    # lclz[i] = mc.thickness_hydrostatic(data['p'], data['T'], data['q']*units('kg/kg'), bottom=sfcp[i]*units.hPa, depth=(sfcp[i]-lclp[i])*units.hPa).magnitude
    
    # pp = np.array([sfcp[i]] + list(p[p<sfcp[i]]))*units.hPa
    # tt = np.array([t2m[i]] + list(data['T'][p<sfcp[i]].magnitude))*units.K
    # dd = np.array([td2m[i]] + list(data['Td'][p<sfcp[i]].magnitude+273.15))*units.K
    # dc = mc.downdraft_cape(pp, tt, dd)
    # dcape[i] = dc[0].magnitude
    # # downp = dc[1].magnitude[0]
    # downT[i] = dc[2].magnitude[0]
    
    
    datap.close()
    datas.close()

#%% Plot bar charts of sounding parameters

figsave = False

wid1 = 0.25
wid2 = 0.2

x0 = np.arange(len(cape))
x1 = [x + wid1 for x in x0]
x2 = [x + wid1 for x in x1]
x3 = [x + wid1 for x in x2]


# Bulk shear
fig,ax = plt.subplots(figsize=(9,5), layout='constrained')

ax.bar(x1, shear01, color='salmon', width=wid2, edgecolor='k', label='0-1 km')
ax.bar(x2, shear03, color='gold', width=wid2, edgecolor='k', label='0-3 km')
ax.bar(x3, shear06, color='skyblue', width=wid2, edgecolor='k', label='0-6 km')

# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('Bulk wind shear (m/s)', fontsize=12)
ax.set_xticks(x2, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], shear01[i], f"{shear01[i]:.0f}", fontsize=12, ha='center', va='bottom')
    ax.text(x2[i], shear03[i], f"{shear03[i]:.0f}", fontsize=12, ha='center', va='bottom')
    ax.text(x3[i], shear06[i], f"{shear06[i]:.0f}", fontsize=12, ha='center', va='bottom')
ax.set_title('Bulk Wind Shear 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')

ax.grid(visible=True, which='major', axis='y', color='lightgray', linewidth=0.75, linestyle='--')
ax.grid(visible=True, which='minor', axis='x', color='lightgray', linewidth=0.75, linestyle='--')
ax.set_axisbelow(True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.set_xlim([0,7])
ax.set_ylim([0,35])
if figsave:
    plt.savefig(fp+'figs/barplot_shear.png', dpi=300)




wid1 = 0.25
wid2 = 0.15

x1 = np.arange(len(cape))
# x2 = [x + wid1 for x in x1]
x1 = [x + 0.5-wid2 for x in x0]
x2 = [x + 0.5+wid2 for x in x0]
# x3 = 
# x4 = [x + wid1 for x in x3]


# SRH
fig,ax = plt.subplots(figsize=(9,5), layout='constrained')

ax.bar(x1, srh01, color='salmon', width=wid1, edgecolor='k', label='0-1 km')
ax.bar(x2, srh03, color='lightskyblue', width=wid1, edgecolor='k', label='0-3 km')


# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('SRH (m2/s2)', fontsize=12)
ax.set_xticks([x + 0.5 for x in x0],
              ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
               '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
               '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], srh01[i], f"{srh01[i]:.0f}", fontsize=12, ha='center', va='bottom')
    ax.text(x2[i], srh03[i], f"{srh03[i]:.0f}", fontsize=12, ha='center', va='bottom')
ax.set_title('Storm-Relative Helicity 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')

ax.grid(visible=True, which='major', axis='y', color='lightgray', linewidth=0.75, linestyle='--')
ax.grid(visible=True, which='minor', axis='x', color='lightgray', linewidth=0.75, linestyle='--')
ax.set_axisbelow(True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(50))
ax.set_xlim([0,7])
ax.set_ylim([0,350])
if figsave:
    plt.savefig(fp+'figs/barplot_srh.png', dpi=300)



# CAPE
wid1 = 0.35
x1 = np.arange(len(cape)) + 0.5

fig,ax = plt.subplots(figsize=(9,5), layout='constrained')

ax.bar(x1, cape, color='lightskyblue', width=wid1, edgecolor='k', label='CAPE')
ax.bar(x1, -1*cin, color='b', width=wid1, edgecolor='k', label='CIN')
# ax.bar(x2, dcape, color='lightskyblue', width=wid1, edgecolor='k', label='DCAPE')
ax.axhline(0, c='k', linewidth=1, linestyle='--')

# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('CAPE/CIN (J/kg)', fontsize=12)
ax.set_xticks(x1, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], cape[i], f"{cape[i]:.0f}", fontsize=12, ha='center', va='bottom')
    ax.text(x1[i], -1*cin[i]-10, f"{-1*cin[i]:.0f}", fontsize=12, ha='center', va='top')
    # ax.text(x2[i], dcape[i], f"{dcape[i]:.0f}", fontsize=10, ha='center', va='bottom')
ax.set_title('CAPE/CIN 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')

ax.grid(visible=True, which='major', axis='y', color='lightgray', linewidth=0.75, linestyle='--')
ax.grid(visible=True, which='minor', axis='x', color='lightgray', linewidth=0.75, linestyle='--')
ax.set_axisbelow(True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(250))
ax.set_xlim([0,7])
ax.set_ylim([-350,2100])
if figsave:
    plt.savefig(fp+'figs/barplot_cape.png', dpi=300)



# LCL pressure and height
wid1 = 0.35
x1 = np.arange(len(cape)) + 0.5
# inverted axis - plot bars=bar(x, anchor-data, ... bottom=data) -> invert_yaxis() -> for b in bars, b.sticky_edges.y[:] = [anchor]

fig,ax = plt.subplots(figsize=(9,5), layout='constrained')

b1 = ax.bar(x1, 1000-lclp, color='tan', width=wid1, edgecolor='gray', bottom=lclp)
bars = ax.bar(x1, sfcp-lclp, color='lightskyblue', width=wid1, edgecolor='k', label='LCL pressure', bottom=lclp)
# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('Pressure (hPa)', fontsize=12)
ax.set_xticks(x1, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], lclp[i], f"{lclp[i]:.0f}", fontsize=12, ha='center', va='bottom')
    ax.text(x1[i], sfcp[i], f"{sfcp[i]:.0f}", fontsize=12, ha='center', va='bottom')
    # ax.plot([x1[i]-wid1, x1[i]+wid1], [sfcp[i], sfcp[i]], 'gray', linewidth=1)
ax.set_title('LCL Pressure 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')

ax.grid(visible=True, which='major', axis='y', color='lightgray', linewidth=0.75, linestyle='--')
ax.grid(visible=True, which='minor', axis='x', color='lightgray', linewidth=0.75, linestyle='--')
ax.set_axisbelow(True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(25))
ax.set_xlim([0,7])
ax.invert_yaxis()
for b in bars:
    b.sticky_edges.y[:] = [sfcp]
for b in b1:
    b.sticky_edges.y[:] = [1000]
ax.set_ylim([1000,800])
if figsave:
    plt.savefig(fp+'figs/barplot_lcl_pressure.png', dpi=300)



fig,ax = plt.subplots(figsize=(9,5), layout='constrained')

bars = ax.bar(x1, lclz, color='lightskyblue', width=wid1, edgecolor='k', label='LCL height')
# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('Height (m)', fontsize=12)
ax.set_xticks(x1, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], lclz[i], f"{lclz[i]:.0f} m", fontsize=12, ha='center', va='bottom')
ax.set_title('LCL Height 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')

ax.grid(visible=True, which='major', axis='y', color='lightgray', linewidth=0.75, linestyle='--')
ax.grid(visible=True, which='minor', axis='x', color='lightgray', linewidth=0.75, linestyle='--')
ax.set_axisbelow(True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(250))
ax.set_xlim([0,7])
ax.set_ylim([0,2000])
if figsave:
    plt.savefig(fp+'figs/barplot_lcl_height.png', dpi=300)



# Lapse rate 700-500 mb
wid1 = 0.35
x1 = np.arange(len(cape)) + 0.5

fig,ax = plt.subplots(figsize=(9,5), layout='constrained')

bars = ax.bar(x1, lr, color='lightskyblue', width=wid1, edgecolor='k', label='Lapse rate')
# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('Lapse rate (K/km)', fontsize=12)
ax.set_xticks(x1, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], lr[i], f"{lr[i]:.1f}", fontsize=12, ha='center', va='bottom')
ax.set_title('700-500 mb Lapse Rate 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')

ax.grid(visible=True, which='major', axis='y', color='lightgray', linewidth=0.75, linestyle='--')
ax.grid(visible=True, which='minor', axis='x', color='lightgray', linewidth=0.75, linestyle='--')
ax.set_axisbelow(True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.set_xlim([0,7])
ax.set_ylim([5,8])
if figsave:
    plt.savefig(fp+'figs/barplot_lapserate.png', dpi=300)



# 2-m dewpoint depression
wid1 = 0.35
x1 = np.arange(len(cape)) + 0.5

fig,ax = plt.subplots(figsize=(9,5), layout='constrained')

bars = ax.bar(x1, t2m-td2m, color='lightskyblue', width=wid1, edgecolor='k', label='T - Td')
# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('T - Td (C)', fontsize=12)
ax.set_xticks(x1, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], t2m[i]-td2m[i], f"{t2m[i]-td2m[i]:.1f} C", fontsize=12, ha='center', va='bottom')
ax.set_title('2-m Dewpoint Depression 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')

ax.grid(visible=True, which='major', axis='y', color='lightgray', linewidth=0.75, linestyle='--')
ax.grid(visible=True, which='minor', axis='x', color='lightgray', linewidth=0.75, linestyle='--')
ax.set_axisbelow(True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(2))
ax.set_xlim([0,7])
ax.set_ylim([0,12])
if figsave:
    plt.savefig(fp+'figs/barplot_dewpt_depression.png', dpi=300)



# DCAPE
wid1 = 0.35
x1 = np.arange(len(cape)) + 0.5

fig,ax = plt.subplots(figsize=(9,5), layout='constrained')

ax.bar(x1, dcape, color='lightskyblue', width=wid1, edgecolor='k', label='DCAPE')
# ax.axhline(0, c='k', linewidth=1, linestyle='--')

# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('DCAPE (J/kg)', fontsize=12)
ax.set_xticks(x1, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], dcape[i], f"{dcape[i]:.0f}", fontsize=12, ha='center', va='bottom')
ax.set_title('Downdraft CAPE 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')

ax.grid(visible=True, which='major', axis='y', color='lightgray', linewidth=0.75, linestyle='--')
ax.grid(visible=True, which='minor', axis='x', color='lightgray', linewidth=0.75, linestyle='--')
ax.set_axisbelow(True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(250))
ax.set_xlim([0,7])
ax.set_ylim([0,1500])
if figsave:
    plt.savefig(fp+'figs/barplot_dcape.png', dpi=300)




# Predicted cold pool temperature deficit
wid1 = 0.35
x1 = np.arange(len(cape)) + 0.5

fig,ax = plt.subplots(figsize=(9,5), layout='constrained')

bars = ax.bar(x1, t2m-downT, color='lightskyblue', width=wid1, edgecolor='k', label='2-m temperature deficit', bottom=downT-t2m)
# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('T(downdraft) - T(env)  (C)', fontsize=12)
ax.set_xticks(x1, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], downT[i]-t2m[i], f"{downT[i]-t2m[i]:.1f} C", fontsize=12, ha='center', va='bottom')
ax.set_title("Predicted Cold Pool Temperature Deficit 1 H Prior to First Tornado/Severe Report", fontsize=14)
ax.legend(fontsize=12, loc='upper left')

ax.grid(visible=True, which='major', axis='y', color='lightgray', linewidth=0.75, linestyle='--')
ax.grid(visible=True, which='minor', axis='x', color='lightgray', linewidth=0.75, linestyle='--')
ax.set_axisbelow(True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(2))
ax.set_xlim([0,7])

ax.invert_yaxis()
for b in bars:
    b.sticky_edges.y[:] = [0]
ax.set_ylim([0,-16])
if figsave:
    plt.savefig(fp+'figs/barplot_coldpool.png', dpi=300)





plt.show()


#%% Box plots of sounding variable distributions along storm path

from era5utils import *

fp = "C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/tor_outbreaks/"

# yyyy = 2025
# mm = 6
# dd = 23

leadtime = 1

dbfile = open(fp+"tornado_locs.pkl", 'rb')
locs = pickle.load(dbfile)
dbfile.close()

events = ["20210811", "20250623", "20220530", "20220521", "20260630", "20250724", "20260703"]

# lats = []
# lons = []
# timt = []
# yyyyt = []
# mmt = []
# ddt = []
# hht = []

data_all = dict()

for i in range(len(events)):
    # tors = locs[events[i]]
    
    # lats = []
    # lons = []
    # timt = []

    # for k in range(len(tors)):
    #     tor = tors[f"loc{k+1}"]
    #     lats.append(tors[f"loc{k+1}"]['lat'])
    #     lons.append(tors[f"loc{k+1}"]['lon'])
    #     [yyyyt,mmt,ddt] = tors[f"loc{k+1}"]['date_ymd']
    #     minute = float(tors[f"loc{k+1}"]['time_utc'][2:])
    #     if minute <= 30:
    #         hht = float(tors[f"loc{k+1}"]['time_utc'][:2]) - leadtime
    #     else:
    #         hht = float(tors[f"loc{k+1}"]['time_utc'][:2]) - leadtime + 1
        
    #     if (hht < 0) or (hht > 23):
    #         yyyyt,mmt,ddt,hht = correct_datetime(yyyyt,mmt,ddt,hht)
        
    #     timt.append(f"{yyyyt}-{mmt:02.0f}-{ddt:02.0f}T{hht:02.0f}:00:00.000000000")
    
    
    [yyyyt,mmt,dayt] = locs[events[i]][f"loc1"]['date_ymd']
    
    lats = latpoints[events[i]]
    lons = lonpoints[events[i]]
    timt = []
    
    for k in range(len(lats)):
        ddt = daypoints[events[i]][k]
        hht = hourpoints[events[i]][k]
        timt.append(f"{yyyyt}-{mmt:02.0f}-{ddt:02.0f}T{hht:02.0f}:00:00.000000000")
    
    
    n = np.arange(len(lats))
    
    t2m = np.zeros((len(lats),))
    td2m = np.zeros((len(lats),))
    cape = np.zeros((len(lats),))
    cin = np.zeros((len(lats),))
    shear01 = np.zeros((len(lats),))
    shear03 = np.zeros((len(lats),))
    shear06 = np.zeros((len(lats),))
    srh01 = np.zeros((len(lats),))
    srh03 = np.zeros((len(lats),))
    lclp = np.zeros((len(lats),))
    lr = np.zeros((len(lats),))
    sfcp = np.zeros((len(lats),))
    lclz = np.zeros((len(lats),))
    dcape = np.zeros((len(lats),))
    downT = np.zeros((len(lats),))
    
    
    for lat,lon,tim,k in zip(lats,lons,timt,n):
        yyyyt = tim[:4]
        mmt = tim[5:7]
        ddt = tim[8:10]
        hht = tim[12:14]
        
        fn_preslev = f"era5_{yyyyt}{mmt}{ddt}_preslevs.nc"
        fn_singlev = f"era5_{yyyyt}{mmt}{ddt}_singlevs.nc"
        
        datap = xr.open_dataset(fp+fn_preslev)
        datas = xr.open_dataset(fp+fn_singlev)
        
        latitude = datap['latitude'][:].values
        longitude = datap['longitude'][:].values
        
        lati = np.argmin(np.abs(latitude-lat))
        loni = np.argmin(np.abs(longitude-lon))
        latt = latitude[lati-1:lati+2]
        lont = longitude[loni-1:loni+2]
        
        # p,z,T,q,theta,Td,u,v,u10,v10,speed,direc,cape,cin,sfcp,orog,q2m,theta2m,td2m,t2m,leftm,meanm,rightm,parcel_prof,lcl_pressure,lcl_temperature = extract_data(latt,lont,timt,dsp,dss)
        data = extract_data(latt,lont,tim,datap,datas)
        
        datap.close()
        datas.close()
        
        p = data['p'].magnitude
        t2m[k] = data['t2m'].magnitude
        td2m[k] = data['td2m'].magnitude
        cape[k] = data['cape']
        cin[k] = data['cin']
        lclp[k] = data['lcl_pressure'].magnitude
        shear06[k] = data['shear06'].magnitude
        shear03[k] = data['shear03'].magnitude
        shear01[k] = data['shear01'].magnitude
        srh03[k] = data['srh03'].magnitude
        srh01[k] = data['srh01'].magnitude
        lr[k] = data['lapse_rate']
        sfcp[k] = data['sfcp']
        lclz[k] = data['lcl_height'].magnitude
        dcape[k] = data['dcape']
        downT[k] = data['downT'].magnitude
        
    dat = dict(p=p, sfcp=sfcp, t2m=t2m, td2m=td2m, cape=cape, cin=cin, lclp=lclp, lclz=lclz, lr=lr,
               shear01=shear01, shear03=shear03, shear06=shear06, srh01=srh01, srh03=srh03, dcape=dcape, downT=downT)
    
    data_all[events[i]] = dat



#%% Make box plots

labels = ['11 Aug 2021\n Outbreak ',
          '23 Jun 2025\n Outbreak ',
          '30 May 2022\n Sub-outbreak ',
          '21 May 2022\n Sub-outbreak ',
          '30 Jun 2026\n Sub-outbreak ',
          '24 Jul 2025\n Null event ',
          '3 July 2026\n Null event ']

labels_none = ['', '', '', '', '', '', '']
bw = 0.25
# colors = ['lightpink', 'salmon', 'crimson']

# colors = 


# fig,ax = plt.subplots(figsize=(8,6), layout='constrained')

# b = ax.boxplot([shear01, shear03, shear06], tick_labels=['0-1 km shear', '0-3 km shear', '0-6 km shear'], patch_artist=True, positions=[0,1,4])

# for patch,c in zip(b['boxes'],colors):
#     patch.set_facecolor(c)

# ax.set_ylabel('Shear [m/s]')
# ax.set_title('Bulk wind difference')
# ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)



shear01_all = [data_all[events[i]]['shear01'] for i in range(len(events))]
shear03_all = [data_all[events[i]]['shear03'] for i in range(len(events))]
shear06_all = [data_all[events[i]]['shear06'] for i in range(len(events))]
cape_all = [data_all[events[i]]['cape'] for i in range(len(events))]
cin_all = [-1*data_all[events[i]]['cin'] for i in range(len(events))]
lclz_all = [data_all[events[i]]['lclz'] for i in range(len(events))]
srh01_all = [data_all[events[i]]['srh01'] for i in range(len(events))]
srh03_all = [data_all[events[i]]['srh03'] for i in range(len(events))]
t2m_all = [data_all[events[i]]['t2m'] for i in range(len(events))]
td2m_all = [data_all[events[i]]['td2m'] for i in range(len(events))]
downT_all = [data_all[events[i]]['downT'] for i in range(len(events))]
tdepr_all = [data_all[events[i]]['t2m'] - data_all[events[i]]['td2m'] for i in range(len(events))]
cpt_all = [data_all[events[i]]['downT'] - data_all[events[i]]['t2m'] for i in range(len(events))]
dcape_all = [data_all[events[i]]['dcape'] for i in range(len(events))]
lr_all = [data_all[events[i]]['lr'] for i in range(len(events))]




fig,ax = plt.subplots(figsize=(9,5), layout='constrained')
b = ax.boxplot(shear01_all, tick_labels=labels, patch_artist=True, positions=np.arange(len(events))/2, widths=bw)
for patch in b['boxes']:
    patch.set_facecolor('lightskyblue')
# b3 = ax.boxplot(shear03_all, tick_labels=labels_none, patch_artist=True, positions=np.arange(len(events))+0.2, widths=bw)
# for patch in b3['boxes']:
#     patch.set_facecolor('blue')
# b6 = ax.boxplot(shear06_all, tick_labels=labels_none, patch_artist=True, positions=np.arange(len(events))+0.4, widths=bw)
# for patch in b6['boxes']:
#     patch.set_facecolor('indigo')
# ax.legend([b1['boxes'][0], b3['boxes'][0], b6['boxes'][0]], ['0-1 km', '0-3 km', '0-6 km'], loc='upper right')
ax.set_ylabel('Shear [m/s]')
ax.set_title('0-1 km bulk wind difference')
ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)


fig,ax = plt.subplots(figsize=(9,5), layout='constrained')
b = ax.boxplot(shear03_all, tick_labels=labels, patch_artist=True, positions=np.arange(len(events))/2, widths=bw)
for patch in b['boxes']:
    patch.set_facecolor('lightskyblue')
ax.set_ylabel('Shear [m/s]')
ax.set_title('0-3 km bulk wind difference')
ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)


fig,ax = plt.subplots(figsize=(9,5), layout='constrained')
b = ax.boxplot(shear06_all, tick_labels=labels, patch_artist=True, positions=np.arange(len(events))/2, widths=bw)
for patch in b['boxes']:
    patch.set_facecolor('lightskyblue')
ax.set_ylabel('Shear [m/s]')
ax.set_title('0-6 km bulk wind difference')
ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)




fig,ax = plt.subplots(figsize=(9,5), layout='constrained')
b = ax.boxplot(cape_all, tick_labels=labels, patch_artist=True, positions=np.arange(len(events))/2, widths=bw)
for patch in b['boxes']:
    patch.set_facecolor('lightskyblue')
ax.set_ylabel('CAPE [J/kg]')
ax.set_title('CAPE')
ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)


fig,ax = plt.subplots(figsize=(9,6), layout='constrained')
b = ax.boxplot(cin_all, tick_labels=labels, patch_artist=True, positions=np.arange(len(events))/2, widths=bw)
for patch in b['boxes']:
    patch.set_facecolor('lightskyblue')
ax.set_ylabel('CIN [J/kg]')
ax.set_title('CIN')
ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)


fig,ax = plt.subplots(figsize=(9,5), layout='constrained')
b = ax.boxplot(dcape_all, tick_labels=labels, patch_artist=True, positions=np.arange(len(events))/2, widths=bw)
for patch in b['boxes']:
    patch.set_facecolor('lightskyblue')
ax.set_ylabel('DCAPE [J/kg]')
ax.set_title('Downdraft CAPE')
ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)


fig,ax = plt.subplots(figsize=(9,5), layout='constrained')
b = ax.boxplot(lr_all, tick_labels=labels, patch_artist=True, positions=np.arange(len(events))/2, widths=bw)
for patch in b['boxes']:
    patch.set_facecolor('lightskyblue')
ax.set_ylabel('Lapse rate [K/km]')
ax.set_title('Mid-level lapse rate (700-500 mb)')
ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)




fig,ax = plt.subplots(figsize=(9,6), layout='constrained')
b = ax.boxplot(srh01_all, tick_labels=labels, patch_artist=True, positions=np.arange(len(events))/2, widths=bw)
for patch in b['boxes']:
    patch.set_facecolor('lightskyblue')
ax.set_ylabel('SRH [m2/s2]')
ax.set_title('0-1 km SRH')
ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)


fig,ax = plt.subplots(figsize=(9,6), layout='constrained')
b = ax.boxplot(srh03_all, tick_labels=labels, patch_artist=True, positions=np.arange(len(events))/2, widths=bw)
for patch in b['boxes']:
    patch.set_facecolor('lightskyblue')
ax.set_ylabel('SRH [m2/s2]')
ax.set_title('0-3 km SRH')
ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)




fig,ax = plt.subplots(figsize=(9,6), layout='constrained')
b = ax.boxplot(lclz_all, tick_labels=labels, patch_artist=True, positions=np.arange(len(events))/2, widths=bw)
for patch in b['boxes']:
    patch.set_facecolor('lightskyblue')
ax.set_ylabel('Height [m]')
ax.set_title('LCL height')
ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)


fig,ax = plt.subplots(figsize=(9,6), layout='constrained')
b = ax.boxplot(tdepr_all, tick_labels=labels, patch_artist=True, positions=np.arange(len(events))/2, widths=bw)
for patch in b['boxes']:
    patch.set_facecolor('lightskyblue')
ax.set_ylabel('Temperature [C]')
ax.set_title('Dewpoint depression')
ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)


fig,ax = plt.subplots(figsize=(9,6), layout='constrained')
b = ax.boxplot(cpt_all, tick_labels=labels, patch_artist=True, positions=np.arange(len(events))/2, widths=bw)
for patch in b['boxes']:
    patch.set_facecolor('lightskyblue')
ax.set_ylabel('Temperature [C]')
ax.set_title('Predicted cold pool temperature deficit')
ax.grid(visible=True, which='both', axis='y', color='lightgray', linewidth=0.75)




plt.show()



#%%

fp = "C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/tor_outbreaks/"

# yyyy = 2025
# mm = 6
# dd = 23

leadtime = 1

dbfile = open(fp+"tornado_locs.pkl", 'rb')
locs = pickle.load(dbfile)
dbfile.close()

events = ["20210811", "20250623", "20220530", "20220521", "20260630", "20250724", "20260703"]


# August 11 2021
lats_1 = [46.5, 46.25, 46.5, 46.25, 46.5, 
          46.4797, 46.5449, 46.5274, 46.4657, 46.25, 46.25, 46.25, 46.5,
          46.3870, 46.6557, 46.25, 46.5, 46.25, 46.25, 46.5, 46.5]
lons_1 = [-84.25, -84.5, -84.0, -84.75, -84.5, 
          -83.6403, -83.5907, -83.2889, -83.2155, -84.0, -83.5, -83.0, -83.0,
          -82.7493, -82.4617, -82.5, -82.5, -81.75, -82.25, -82.25, -82.0]
hours_1 = [19, 19, 19, 19, 19,
           20, 20, 20, 20, 20, 20, 20, 20,
           21, 21, 21, 21, 21, 21, 21, 21]
days_1 = [11, 11, 11, 11, 11,
          11, 11, 11, 11, 11, 11, 11, 11,
          11, 11, 11, 11, 11, 11, 11, 11]
latlont1 = []
for lat,lon,t in zip(lats_1,lons_1,hours_1):
    latr = np.round(lat/0.25) * 0.25
    lonr = np.round(lon/0.25) * 0.25
    latlont1.append([t,latr,lonr])


# June 23 2025
lats_2 = [48.3975, 48.4647, 48.0, 48.25, 48.5, 48.5, 48.75,
          48.5, 48.25, 48.25, 48.5, 48.25, 48.0, 48.0, 
          48.2033, 47.9149, 48.0, 48.0, 48.25, 48.0, 48.25, 47.75,
          47.5, 47.5, 47.5, 47.5, 47.75, 47.75, 47.75,
          47.25, 47.25, 47.25, 47.0, 47.0, 47.0, 47.0, 
          46.7872, 46.7757, 46.75, 46.75, 47.0, 46.5, 46.5]
lons_2 = [-75.5880, -75.4900, -76.75, -76.75, -76.75, -75.25, -75.25,
          -75.0, -74.5, -74.75, -74.25, -74.25, -74.75, -74.5, 
          -73.5234, -73.5056, -73.25, -73.0, -73.25, -73.75, -73.75, -73.25,
          -72.0, -72.25, -72.75, -72.5, -72.75, -72.5, -72.25,
          -72.0, -71.75, -71.5, -72.25, -71.75, -71.5, -71.25, 
          -70.7827, -70.4072, -71.25, -71.0, -70.75, -70.75, -70.5]
hours_2 = [20, 20, 20, 20, 20, 20, 20,
           21, 21, 21, 21, 21, 21, 21,
           22, 22, 22, 22, 22, 22, 22, 22,
           23, 23, 23, 23, 23, 23, 23,
           0, 0, 0, 0, 0, 0, 0,
           1, 1, 1, 1, 1, 1, 1]
days_2 = [23, 23, 23, 23, 23, 23, 23,
          23, 23, 23, 23, 23, 23, 23,
          23, 23, 23, 23, 23, 23, 23, 23,
          23, 23, 23, 23, 23, 23, 23,
          24, 24, 24, 24, 24, 24, 24,
          24, 24, 24, 24, 24, 24, 24]
latlont2 = []
for lat,lon,t in zip(lats_2,lons_2,hours_2):
    latr = np.round(lat/0.25) * 0.25
    lonr = np.round(lon/0.25) * 0.25
    latlont2.append([t,latr,lonr])


# May 30 2022
lats_3 = [48.6148, 48.6699, 48.9114, 48.25, 48.5, 48.5, 48.75, 
          49.3689, 48.6757]
lons_3 = [-93.5304, -93.5219, -93.5778, -93.5, -94.0, -93.75, -93.25, 
          -93.0754, -92.2308]
hours_3 = [0, 0, 0, 0, 0, 0, 0,
           1, 1]
days_3 = [31, 31, 31, 31, 31, 31, 31,
          31, 31]
latlont3 = []
for lat,lon,t in zip(lats_3,lons_3,hours_3):
    latr = np.round(lat/0.25) * 0.25
    lonr = np.round(lon/0.25) * 0.25
    latlont3.append([t,latr,lonr])


# May 21 2022
lats_4 = [42.5, 42.25, 42.5,
          42.5, 42.75, 42.75, 43.0, 42.75, 43.0,
          43.0179, 42.9217, 42.75, 43.0, 43.25, 43.0, 43.25, 43.25,
          43.5, 43.25, 43.75, 43.75, 43.5, 43.75,
          44.1058, 44.1755, 44.0, 43.75, 44.0, 44.0, 44.25]
lons_4 = [-83.0, -83.75, -83.25,
          -83.0, -82.75, -82.5, -82.25, -82.25, -82.0,
          -81.2216, -81.1977, -81.75, -81.5, -81.25, -81.0, -80.75, -81.0,
          -80.5, -80.5, -80.5, -80.0, -80.25, -80.25,
          -79.1458, -78.7722, -79.5, -79.5, -79.0, -78.5, -78.5]
hours_4 = [13, 13, 13,
           14, 14, 14, 14, 14, 14,
           15, 15, 15, 15, 15, 15, 15, 15,
           16, 16, 16, 16, 16, 16,
           17, 17, 17, 17, 17, 17, 17]
days_4 = [21, 21, 21,
          21, 21, 21, 21, 21, 21,
          21, 21, 21, 21, 21, 21, 21, 21,
          21, 21, 21, 21, 21, 21,
          21, 21, 21, 21, 21, 21, 21]
latlont4 = []
for lat,lon,t in zip(lats_4,lons_4,hours_4):
    latr = np.round(lat/0.25) * 0.25
    lonr = np.round(lon/0.25) * 0.25
    latlont4.append([t,latr,lonr])


# June 30 2026
lats_5 = [45.6315, 44.6890,
          44.4640, 44.3021, 44.3303,
          44.0, 43.75]
lons_5 = [-77.8432, -76.9664,
          -76.7289, -76.5386, -76.9268,
          -75.5, -75.5]
hours_5 = [15, 15,
           16, 16, 16,
           17, 17]
days_5 = [30, 30,
          30, 30, 30,
          30, 30]
latlont5 = []
for lat,lon,t in zip(lats_5,lons_5,hours_5):
    latr = np.round(lat/0.25) * 0.25
    lonr = np.round(lon/0.25) * 0.25
    latlont5.append([t,latr,lonr])


# July 24 2025
lats_6 = [43.47, 44.25, 45.00]
lons_6 = [-81.18, -80.50, -79.25]
hours_6 = [23, 23, 23]
days_6 = [24, 24, 24]
latlont6 = []
for lat,lon,t in zip(lats_6,lons_6,hours_6):
    latr = np.round(lat/0.25) * 0.25
    lonr = np.round(lon/0.25) * 0.25
    latlont6.append([t,latr,lonr])


# July 3 2026
lats_7 = [42.2766,
          42.2393,
          42.4226,
          42.7505]
lons_7 = [-83.7341,
          -83.0576,
          -82.1339,
          -81.7004]
hours_7 = [21,
           22,
           23,
           0]
days_7 = [3,
          3,
          3,
          4]
latlont7 = []
for lat,lon,t in zip(lats_7,lons_7,hours_7):
    latr = np.round(lat/0.25) * 0.25
    lonr = np.round(lon/0.25) * 0.25
    latlont7.append([t,latr,lonr])




latpoints = {"20210811":lats_1, "20250623":lats_2, "20220530":lats_3, "20220521":lats_4, "20260630":lats_5, "20250724":lats_6, "20260703":lats_7}
lonpoints = {"20210811":lons_1, "20250623":lons_2, "20220530":lons_3, "20220521":lons_4, "20260630":lons_5, "20250724":lons_6, "20260703":lons_7}
hourpoints = {"20210811":hours_1, "20250623":hours_2, "20220530":hours_3, "20220521":hours_4, "20260630":hours_5, "20250724":hours_6, "20260703":hours_7}
daypoints = {"20210811":days_1, "20250623":days_2, "20220530":days_3, "20220521":days_4, "20260630":days_5, "20250724":days_6, "20260703":days_7}























