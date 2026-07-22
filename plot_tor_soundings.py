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



#%% 

capes1 = {'0h':1111.8, '1h':1076.5, '2h':544.8}
cins1 = {'0h':-167.3, '1h':-52.1, '2h':-414.3}
shear1s1 = {'0h':12.6, '1h':10.3, '2h':10.1}
shear3s1 = {'0h':20.1, '1h':16.5, '2h':15.6}
shear6s1 = {'0h':26.3, '1h':22.4, '2h':22.0}
srh3s1 = {'0h':281.2, '1h':198.5, '2h':150.2}
srh1s1 = {'0h':162.5, '1h':111.9, '2h':69.5}
lcls1 = {'0h':933.3, '1h':904.5, '2h':877.8}


capes2 = {'0h':1177.3, '1h':1131.4, '2h':1158.6}
cins2 = {'0h':-73.6, '1h':-114.9, '2h':-108.0}
shear1s2 = {'0h':12.8, '1h':12.4, '2h':11.6}
shear3s2 = {'0h':14.2, '1h':14.9, '2h':15.1}
shear6s2 = {'0h':15.9, '1h':17.1, '2h':16.5}
srh3s2 = {'0h':165.1, '1h':221.3, '2h':208.6}
srh1s2 = {'0h':136.7, '1h':141.4, '2h':117.3}
lcls2 = {'0h':901.8, '1h':889.4, '2h':871.0}


capes3 = {'0h':725.2, '1h':762.6, '2h':904.6}
cins3 = {'0h':-51.0, '1h':-49.1, '2h':-54.6}
shear1s3 = {'0h':14.1, '1h':13.9, '2h':14.6}
shear3s3 = {'0h':18.6, '1h':19.4, '2h':21.5}
shear6s3 = {'0h':32.4, '1h':30.8, '2h':29.4}
srh3s3 = {'0h':359.6, '1h':362.2, '2h':317.5}
srh1s3 = {'0h':243.8, '1h':219.2, '2h':182.8}
lcls3 = {'0h':879.1, '1h':879.3, '2h':903.5}


capes4 = {'0h':1044.0, '1h':1286.7, '2h':1149.3}
cins4 = {'0h':-125.3, '1h':-121.1, '2h':-183.0}
shear1s4 = {'0h':10.5, '1h':9.7, '2h':9.1}
shear3s4 = {'0h':17.4, '1h':15.9, '2h':16.9}
shear6s4 = {'0h':25.1, '1h':23.2, '2h':24.2}
srh3s4 = {'0h':110.1, '1h':114.8, '2h':123.2}
srh1s4 = {'0h':52.5, '1h':76.3, '2h':78.7}
lcls4 = {'0h':915.3, '1h':903.9, '2h':903.0}


capes5 = {'0h':1561.8, '1h':1473.4, '2h':1238.1}
cins5 = {'0h':-103.2, '1h':-263.0, '2h':-304.0}
shear1s5 = {'0h':11.6, '1h':9.8, '2h':8.5}
shear3s5 = {'0h':13.3, '1h':11.9, '2h':13.1}
shear6s5 = {'0h':17.2, '1h':18.5, '2h':17.5}
srh3s5 = {'0h':92.3, '1h':112.6, '2h':117.6}
srh1s5 = {'0h':73.2, '1h':89.9, '2h':82.6}
lcls5 = {'0h':908.3, '1h':870.3, '2h':842.1}


capes6 = {'-1h':505.3, '0h':530.6, '1h':995.1, '2h':1432.5}
cins6 = {'-1h':-74.7, '0h':-31.2, '1h':-63.8, '2h':-75.0}
shear1s6 = {'-1h':20.3, '0h':22.3, '1h':17.4, '2h':13.4}
shear3s6 = {'-1h':20.3, '0h':21.8, '1h':13.8, '2h':8.1}
shear6s6 = {'-1h':25.1, '0h':26.5, '1h':24.5, '2h':18.9}
srh3s6 = {'-1h':337.5, '0h':385.7, '1h':175.5, '2h':12.6}
srh1s6 = {'-1h':275.5, '0h':325.5, '1h':183.6, '2h':88.4}
lcls6 = {'-1h':938.6, '0h':942.0, '1h':941.4, '2h':935.5}



# April 15 with Michigan tors
capes6_2 = {'-1h':482.1, '0h':497.2, '1h':959.8, '2h':1280.4}
cins6_2 = {'-1h':-62.3, '0h':-34.1, '1h':-80.6, '2h':-118.6}
shear1s6_2 = {'-1h':20.6, '0h':22.5, '1h':19.0, '2h':14.7}
shear3s6_2 = {'-1h':20.5, '0h':22.0, '1h':15.8, '2h':10.1}
shears6_2 = {'-1h':24.7, '0h':26.6, '1h':25.3, '2h':19.8}
srh3s6_2 = {'-1h':364.6, '0h':448.7, '1h':244.1, '2h':30.9}
srh1s6_2 = {'-1h':286.3, '0h':352.5, '1h':224.3, '2h':94.4}
lcls6_2 = {'-1h':938.6, '0h':942.3, '1h':939.7, '2h':935.6}



#%% First tornado only

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
    
    latt = latitude[lati-4:lati+4]
    lont = longitude[loni-4:loni+4]
    
    # p,z,T,q,theta,Td,u,v,u10,v10,speed,direc,cape,cin,sfcp,orog,q2m,theta2m,td2m,t2m,leftm,meanm,rightm,parcel_prof,lcl_pressure,lcl_temperature = extract_data(latt,lont,timt,dsp,dss)
    data = extract_data(latt, lont, timt, datap, datas)
    
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

#%%

figsave = True

wid1 = 0.25
wid2 = 0.2

x1 = np.arange(len(cape))
x2 = [x + wid1 for x in x1]
x3 = [x + wid1 for x in x2]


# Bulk shear
fig,ax = plt.subplots(figsize=(8,5), layout='constrained')

ax.bar(x1, shear01, color='salmon', width=wid2, edgecolor='k', label='0-1 km')
ax.bar(x2, shear03, color='gold', width=wid2, edgecolor='k', label='0-3 km')
ax.bar(x3, shear06, color='skyblue', width=wid2, edgecolor='k', label='0-6 km')

# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('Bulk wind shear (m/s)', fontsize=12)
ax.set_xticks(x2, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], shear01[i], f"{shear01[i]:.0f}", fontsize=10, ha='center', va='bottom')
    ax.text(x2[i], shear03[i], f"{shear03[i]:.0f}", fontsize=10, ha='center', va='bottom')
    ax.text(x3[i], shear06[i], f"{shear06[i]:.0f}", fontsize=10, ha='center', va='bottom')
ax.set_title('Bulk Wind Shear 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper right')
if figsave:
    plt.savefig(fp+'figs/barplot_shear.png', dpi=300)



wid1 = 0.2
wid2 = 0.15

x1 = np.arange(len(cape))*1.1
# x2 = [x + wid1 for x in x1]
x2 = [x + wid1+0.05 for x in x1]
# x4 = [x + wid1 for x in x3]


# SRH
fig,ax = plt.subplots(figsize=(8,5), layout='constrained')

ax.bar(x1, srh01, color='salmon', width=wid1, edgecolor='k', label='0-1 km')
ax.bar(x2, srh03, color='lightskyblue', width=wid1, edgecolor='k', label='0-3 km')


# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('SRH (m2/s2)', fontsize=12)
ax.set_xticks((x1+x2)/2,
              ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
               '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
               '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], srh01[i], f"{srh01[i]:.0f}", fontsize=10, ha='center', va='bottom')
    ax.text(x2[i], srh03[i], f"{srh03[i]:.0f}", fontsize=10, ha='center', va='bottom')
ax.set_title('Storm-Relative Helicity 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=11, loc='upper left')
if figsave:
    plt.savefig(fp+'figs/barplot_srh.png', dpi=300)



# CAPE
wid1 = 0.25
x1 = np.arange(len(cape))
x2 = [x + wid1+0.05 for x in x1]

fig,ax = plt.subplots(figsize=(8,5), layout='constrained')

ax.bar(x1, cape, color='lightskyblue', width=wid1, edgecolor='k', label='CAPE')
ax.bar(x1, -1*cin, color='b', width=wid1, edgecolor='k', label='CIN')
ax.axhline(0, c='k', linewidth=1, linestyle='--')

# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('CAPE/CIN (J/kg)', fontsize=12)
ax.set_xticks(x1, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], cape[i], f"{cape[i]:.0f}", fontsize=10, ha='center', va='bottom')
    ax.text(x1[i], -1*cin[i], f"{-1*cin[i]:.0f}", fontsize=10, ha='center', va='top')
ax.set_title('CAPE/CIN 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')
if figsave:
    plt.savefig(fp+'figs/barplot_cape.png', dpi=300)



# LCL
wid1 = 0.25
x1 = np.arange(len(cape))
x2 = [x + wid1+0.05 for x in x1]

fig,ax = plt.subplots(figsize=(8,5), layout='constrained')

bars = ax.bar(x1, 1000-lclp, color='lightskyblue', width=wid1, edgecolor='k', label='LCL pressure', bottom=lclp)
# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('Pressure (hPa)', fontsize=12)
ax.set_xticks(x1, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], lclp[i], f"{lclp[i]:.0f}", fontsize=10, ha='center', va='bottom')
ax.set_title('LCL Pressure 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')

ax.invert_yaxis()
for b in bars:
    b.sticky_edges.y[:] = [1000]
ax.set_ylim([1000,800])
if figsave:
    plt.savefig(fp+'figs/barplot_lcl.png', dpi=300)



# Lapse rate 700-500 mb
wid1 = 0.25
x1 = np.arange(len(cape))
x2 = [x + wid1+0.05 for x in x1]

fig,ax = plt.subplots(figsize=(8,5), layout='constrained')

bars = ax.bar(x1, lr, color='lightskyblue', width=wid1, edgecolor='k', label='Lapse rate')
# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('Lapse rate (K/km)', fontsize=12)
ax.set_xticks(x1, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], lr[i], f"{lr[i]:.1f}", fontsize=10, ha='center', va='bottom')
ax.set_title('700-500 mb lapse rate 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')
ax.set_ylim([5,8])
if figsave:
    plt.savefig(fp+'figs/barplot_lapserate.png', dpi=300)



# 2-m dewpoint depression
wid1 = 0.25
x1 = np.arange(len(cape))
x2 = [x + wid1+0.05 for x in x1]

fig,ax = plt.subplots(figsize=(8,5), layout='constrained')

bars = ax.bar(x1, t2m-td2m, color='lightskyblue', width=wid1, edgecolor='k', label='T - Td')
# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('T - Td (C)', fontsize=12)
ax.set_xticks(x1, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ', '30 Jun 2026\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '3 Jul 2026\n Null event '], fontsize=10)
for i in range(len(x1)):
    ax.text(x1[i], t2m[i]-td2m[i], f"{t2m[i]-td2m[i]:.1f}", fontsize=10, ha='center', va='bottom')
ax.set_title('2-m Dewpoint Depression 1 H Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')
if figsave:
    plt.savefig(fp+'figs/barplot_dewpt_depression.png', dpi=300)


plt.show()


#%%


cape_1 = [1469.8, 998.0, 1096.9, 1061.2, 1583.5, 880.6]
cape_0 = [1184.5, 958.8, 671.8, 840.9, 1526.0, 476.9]

cin_1 = [3.9, 40.7, 38.0, 121.4, 178.6, 70.2]
cin_0 = [115.5, 6.2, 87.4, 115.4, 13.4, 23.5]

shear1_1 = [11.2, 11.2, 14.9, 13.2, 8.3, 18.7]
shear1_0 = [13.9, 10.7, 13.7, 13.6, 9.4, 22.4]

shear3_1 = [16.6, 16.0, 20.7, 16.3, 9.2, 15.1]
shear3_0 = [21.8, 15.7, 21.2, 16.6, 11.7, 22.1]

shear6_1 = [23.7, 16.6, 31.2, 24.0, 15.5, 25.1]
shear6_0 = [29.5, 16.2, 36.8, 24.8, 15.8, 26.5]

srh1_1 = [144.2, 137.4, 177.1, 107.8, 82.8, 210.6]
srh1_0 = [227.0, 142.0, 238.7, 66.1, 13.9, 347.8]

srh3_1 = [216.5, 235.7, 380.6, 124.4, 100.9, 215.5]
srh3_0 = [324.6, 194.6, 404.8, 90.7, 33.3, 422.6]

lcl_1 = [928.5, 878.7, 912.4, 925.5, 841.6, 940.2]
lcl_0 = [946.5, 890.2, 873.7, 928.1, 897.9, 941.6]



wid1 = 0.25
wid2 = 0.2

x1 = np.arange(len(cape_1))
x2 = [x + wid1 for x in x1]
x3 = [x + wid1 for x in x2]


# Bulk shear
fig,ax = plt.subplots(figsize=(8,5), layout='constrained')

ax.bar(x1, shear1_0, color='r', width=wid2, edgecolor='k', label='0-1 km')
ax.bar(x2, shear3_0, color='gold', width=wid2, edgecolor='k', label='0-3 km')
ax.bar(x3, shear6_0, color='b', width=wid2, edgecolor='k', label='0-6 km')

# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('Bulk wind shear (m/s)', fontsize=12)
ax.set_xticks(x2, ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
                   '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ',
                   '24 Jul 2025\n Null event ', '15 Apr 2026\n Null event '], fontsize=11)
ax.set_title('Bulk Wind Shear Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper left')




wid1 = 0.2
wid2 = 0.15

x1 = np.arange(len(cape_1))*1.1
x2 = [x + wid1 for x in x1]
x3 = [x + wid1+0.05 for x in x2]
x4 = [x + wid1 for x in x3]


# SRH
fig,ax = plt.subplots(figsize=(8,5), layout='constrained')

ax.bar(x1, srh1_1, color='pink', width=wid1, edgecolor='k', label='0-1 km (T-1H)')
ax.bar(x2, srh3_1, color='lightskyblue', width=wid1, edgecolor='k', label='0-3 km (T-1H)')
ax.bar(x3, srh1_0, color='r', width=wid1, edgecolor='k', label='0-1 km (T-0H)')
ax.bar(x4, srh3_0, color='b', width=wid1, edgecolor='k', label='0-3 km (T-0H)')

ax.plot([x1, x3], [srh1_1, srh1_0], '--k')
ax.plot([x2, x4], [srh3_1, srh3_0], '--k')

# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('SRH (m2/s2)', fontsize=12)
ax.set_xticks([x + wid2/2+0.05 for x in x2], 
              ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
               '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ',
               '24 Jul 2025\n Null event ', '15 Apr 2026\n Null event '], fontsize=11)
ax.set_title('Storm-Relative Helicity Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=11, loc='upper left')



# CAPE
wid1 = 0.25
x1 = np.arange(len(cape_1))
x2 = [x + wid1+0.05 for x in x1]

fig,ax = plt.subplots(figsize=(8,5), layout='constrained')

ax.bar(x1, cape_1, color='lightskyblue', width=wid1, edgecolor='k', label='T-1H')
ax.bar(x2, cape_0, color='b', width=wid1, edgecolor='k', label='T-0H')

# ax.set_xlabel('Event date', fontsize=12)
ax.set_ylabel('CAPE (J/kg)', fontsize=12)
ax.set_xticks([x + wid1/2 for x in x1], 
              ['11 Aug 2021\n Outbreak ', '23 Jun 2025\n Outbreak ',
               '30 May 2022\n Sub-outbreak ', '21 May 2022\n Sub-outbreak ',
               '24 Jul 2025\n Null event ', '15 Apr 2026\n Null event '], fontsize=11)
ax.set_title('CAPE Prior to First Tornado/Severe Report', fontsize=16)
ax.legend(fontsize=12, loc='upper right')



plt.show()








