# -*- coding: utf-8 -*-
"""
Created on Fri May  1 12:00:28 2026

@author: mschne28
"""

### Import packages ###

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import netCDF4 as nc
import xarray as xr
import pandas as pd
import os
from matplotlib.ticker import MultipleLocator
import pickle

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, BoundaryNorm
from matplotlib import  cm
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors

import metpy.calc as mc
from metpy.cbook import get_test_data
from metpy.plots import Hodograph, SkewT
from metpy.units import units
from metpy.calc import dewpoint_from_relative_humidity,specific_humidity_from_dewpoint
from metpy.calc import potential_temperature,temperature_from_potential_temperature
from metpy.calc import wind_speed, wind_direction,bunkers_storm_motion

import pyart
import cmocean



# import geopandas as gpd
# from shapely.geometry import LineString, Point

#%% Function definitions

# Extract data at location/time
def extract_data(latt,lont,timt,data,datas):
    data01 = data.sel(latitude=latt, longitude=lont, valid_time=timt)
    data02 = datas.sel(latitude=latt, longitude=lont, valid_time=timt)

    p = data01.pressure_level.values * units.hPa
    z = data01['z'].values/9.81
    T = data01['t'].values*units.K
    q = data01['q'].values 
    theta = mc.potential_temperature(p, T)
    Td = mc.dewpoint_from_relative_humidity(T, data01['r']*units.percent)
    u = data01['u'].values*units('m/s')
    v = data01['v'].values*units('m/s')
    speed = mc.wind_speed(u,v)
    direc = mc.wind_direction(u,v)

    cape = data02['cape'].values
    cin  = data02['cin'].values
    sfcp = data02['sp'].values
    orog = data02['z'].values/9.81
    q2m  = mc.specific_humidity_from_dewpoint(sfcp*units.Pa, data02['d2m'].values*units.K)
    theta2m = mc.potential_temperature(sfcp*units.Pa, data02['t2m'].values*units.K)
    td2m = data02['d2m'].values * units.K
    t2m  = data02['t2m'].values * units.K
    u10 = data02['u10'].values * units('m/s')
    v10 = data02['v10'].values * units('m/s')
    
    # Calculate Bunkers Storm Motion:
    [rightm,leftm,meanm] = mc.bunkers_storm_motion(p[z>orog], u[z>orog], v[z>orog], z[z>orog]*units.m)
    
    # Calculate the LCL
    lcl_pressure,lcl_temperature = mc.lcl(p[z>orog][0], T[z>orog][0], Td[z>orog][0])
    
    # Calculate the parcel profile.
    # calclow = 0
    # for i in range(0,len(p)):
    #     if calclow == 0:
    #         if z[i] > orog:
    #             parcel_prof = mc.parcel_profile(p[i:], T[i], Td[i]).to('degC')
    #             #parcel_prof = mc.parcel_profile(p, T[0], Td[0]).to('degC')
    #             calclow = 1
    #             parcelstart = i
    parcel_prof = np.zeros((len(p),))
    parcel_prof[z>orog] = mc.parcel_profile(p[z>orog], T[z>orog][0], Td[z>orog][0]).to('degC')
    parcel_prof[z<=orog] = np.nan
    # parcel_prof = mc.parcel_profile(p, T[0], Td[0]).to('degC')
    
    data_out = {'p':p, 'z':z, 'T':T, 'q':q, 'theta':theta, 'Td':Td, 'u':u, 'v':v, 'u10':u10, 'v10':v10, 'speed':speed, 'direc':direc,
                'cape':cape, 'cin':cin, 'sfcp':sfcp, 'orog':orog, 'q2m':q2m, 'theta2m':theta2m, 'td2m':td2m, 't2m':t2m,
                'leftm':leftm, 'meanm':meanm, 'rightm':rightm, 'parcel_prof':parcel_prof*units.degC, 'lcl_pressure':lcl_pressure, 'lcl_temperature':lcl_temperature}
    return data_out
    # return p,z,T,q,theta,Td,u,v,u10,v10,speed,direc,cape,cin,sfcp,orog,q2m,theta2m,td2m,t2m,leftm,meanm,rightm,parcel_prof,lcl_pressure,lcl_temperature



# Plot skew-T with hodograph -- Mine
def plot_skewT(timt,latt,lont,data,tloc,titlestr=None,figname=None,figfolder=None,figsave=False):
    z = data['z']
    orog = data['orog']
    T = data['T']
    Td = data['Td']
    p = data['p']
    u = data['u']
    v = data['v']
    u10 = data['u10']
    v10 = data['v10']
    leftmover = data['leftm']
    meanwind = data['meanm']
    rightmover = data['rightm']
    lcl_pressure = data['lcl_pressure']
    lcl_temperature = data['lcl_temperature']
    parcel_prof = data['parcel_prof']
    
    name = tloc["name"]
    yyyyt = tloc["date_ymd"][0]
    mmt = tloc["date_ymd"][1]
    ddt = tloc["date_ymd"][2]
    timestr = tloc["time_utc"]
    lattstart = tloc["lat"]
    lontstart = tloc["lon"]
    
    
    
    fig = plt.figure(figsize=(9,9))
    skew = SkewT(fig)

    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-50, 40)

    # Plot special lines
    skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black') # Plot LCL as black dot
    skew.ax.axvline(0, color='c', linestyle='--', linewidth=2) # Plot a zero degree isotherm
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()

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
    if titlestr is not None:
        skew.ax.set_title(titlestr, fontsize=12)
    else:
        skew.ax.set_title(f"{name} , {yyyyt}-{mmt:02.0f}-{ddt:02.0f} T{timestr[:2]}:{timestr[2:]} UTC ({lattstart:.2f}, {lontstart:.2f}) \n ERA5 sounding: {timt[:16]} UTC ({latt:.2f}, {lont:.2f})", fontsize=12)

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

    if figsave is not None:
        if figname is not None:
            plt.savefig(figfolder+figname+'.png', facecolor='white', transparent=False, dpi=300)
        else:
            plt.savefig(figfolder+f"input_sounding_{timt[:13]}_{latt:.2f}_{lont:.2f}.png", facecolor='white', transparent=False, dpi=300)
        # plt.close()
    
    plt.show()
    return fig






# Plot skew-T with hodograph -- Lisa's
def plot_with_hodograph(timt,lont,latt,z,orog,T,Td,p,u,v,u10,v10,leftmover,meanflow,rightmover,lcl_pressure,lcl_temperature,parcel_prof,pngfolder=None):
    # Create a new figure. The dimensions here give a good aspect ratio
    fig = plt.figure(figsize=(9,9))
    skew = SkewT(fig, rotation=30)
    
    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
 #   skew.plot(p[z>orog], temperature_from_potential_temperature(p[z>orog],thetainv[z>orog]), 'b-', linewidth=2)
    skew.plot(p[z>orog], T[z>orog], 'r')
    skew.plot(p[z>orog], Td[z>orog], 'g')
    skew.plot_barbs(p[z>orog], u[z>orog], v[z>orog], xloc=1.1, plot_units=units('kts'))
    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-50, 40)
    
    # Plot LCL as black dot
    skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')
    
    # Plot the parcel profile as a black line
    skew.plot(p[z>orog], parcel_prof, 'k', linewidth=2)
    
    # Shade areas of CAPE and CIN
    skew.shade_cin(p[z>orog], T[z>orog], parcel_prof, Td[z>orog])
    skew.shade_cape(p[z>orog], T[z>orog], parcel_prof)
    
    # Plot a zero degree isotherm
    skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)
    
    # Add the relevant special lines
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()

    plt.title("%s UTC , (lon,lat)=(%.2f,%.2f)" %(timt[:13],lont,latt))
#    plt.imshow([[100,200,300,400,500,600,700,800,900,1000]],cmap=plt.get_cmap('plasma_r'))
    
    # Create a hodograph
    # Create an inset axes object that is 40% width and height of the
    # figure and put it in the upper right hand corner.
    ax_hod = inset_axes(skew.ax, '40%', '40%', loc=1)
    h = Hodograph(ax_hod, component_range=35.)
    h.add_grid(increment=10)
    
    # Get the tab20 colormap from matplotlib
    cmap = plt.get_cmap('tab20')
    # Extract the first 7 colors from it
    colors = [cmap(i) for i in range(7)]
#    colors = sns.cubehelix_palette(5)
#    colors=sns.color_palette("tab20",7)
#    colors=sns.color_palette("rocket",7)
    #colors = ["#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494"]
    my_cmap = ListedColormap(colors, name="my_cmap")
    my_cmap_r = my_cmap.reversed()
    
#    h.plot_colormapped(u[:-1], v[:-1], speed[:-1])  # Plot a line colored by wind speed
    unew = np.zeros((len(u[z>orog])+1))
    vnew = np.zeros((len(v[z>orog])+1))
    znew = np.zeros((len(z[z>orog])+1))
    for k in range(0,len(unew)):
        if k == 0:
            unew[k] = u10.magnitude
            vnew[k] = v10.magnitude
            znew[k] = 10
        else:
            unew[k] = u[z>orog][k-1].magnitude
            vnew[k] = v[z>orog][k-1].magnitude
            znew[k] = z[z>orog][k-1] - orog
    
    hplot = h.plot_colormapped(unew*units("m/s"), vnew*units("m/s"), znew/1000, intervals=[0,1,3,5,7.5,10,12.5], colors=colors)
#    hplot=h.plot_colormapped([u10,u[z>orog]], [v10,v[z>orog]], [10,z[z>orog]-orog],
#                             intervals=[1000,3000,5000,7500,10000],colors=colors)
#    hplot=h.plot_colormapped(u[:-1], v[:-1], p[:-1],intervals=[200,300,500,700,850,1000],colors=colors)
                             #cmap=plt.get_cmap('plasma_r'),
 #                      norm=colors.Normalize(100,1000))  # Plot a line colored by wind speed
    
    cax = ax_hod.inset_axes([0.06, 0.1, 0.88, 0.04])    
    fig.colorbar(hplot, cax=cax, orientation='horizontal',ticks=[0,1,3,5,7.5,10,12.5])#[200, 300,500, 700,850,1000])

    ax_hod.scatter(leftmover[0].magnitude, leftmover[1].magnitude, s=10, color='b')
    ax_hod.scatter(rightmover[0].magnitude, rightmover[1].magnitude, s=10, color='b')
    ax_hod.scatter(meanflow[0].magnitude, meanflow[1].magnitude, s=10, color='b')
    ax_hod.scatter(u10.magnitude, v10.magnitude, s=25, color='k')

    if pngfolder is not None:
        plt.savefig("%s/input_sounding_%s_%s_%s.png"  %(pngfolder,timt[:13],lont,latt), facecolor='white', transparent=False)
        plt.close()
    else:
        plt.show()

    # Show the plot
    # plt.show()
    return fig



def InterpolateToHeightAboveGround(z, orog, param, height):
    nlevs = z.shape[0] 
    nx = z.shape[1]
    ny = z.shape[2]
    
    z2 = np.zeros((int(nlevs),int(nx),int(ny)))
    param_int = np.zeros((int(nx),int(ny)))
    for k in range(0,nlevs):
        z2[k,:,:] = z[k,:,:]-orog
    
    for i in range(0,nx):
        for j in range(0,ny):
            z2idx = np.argmax(np.where(z2[:,i,j] < height))
            param_below = param[int(z2idx),i,j]
            param_above = param[int(z2idx+1),i,j]
            height_below = z2[int(z2idx),i,j]
            height_above = z2[int(z2idx+1),i,j]
            param_int[i,j] = (param_above-param_below) / (height_above-height_below) * (height-height_below) + param_below

    return param_int



# Plot maps
def PlotMoistureMap1(data,datas,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,case,n,pngfolder=None):
#    data01 = data.sel(pressure_level=925, valid_time=timt)
    data01 = data.sel(pressure_level=850, valid_time=timt)
    data03 = data.sel(pressure_level=500, valid_time=timt)
    data02 = datas.sel(valid_time=timt)
    data00 = data.sel(valid_time=timt)
    orog = data02['z'].values/9.81
    z00 = data00['z'].values/9.81
    u00 = data00['u'].values
    v00 = data00['v'].values
    z00 = data00['z'].values/9.81
    q00 = data00['q'].values
    
    u500 = InterpolateToHeightAboveGround(z00, orog, np.array(u00), 500)
    v500 = InterpolateToHeightAboveGround(z00, orog, np.array(v00), 500)
    q500 = InterpolateToHeightAboveGround(z00, orog, np.array(q00), 500)

#    q = data01["q"] 
#    z = data01["z"] / (9.81*10) #in gpdm
#    uu = data01["u"]
#    vv = data01["v"]
    t500 = data03['t'].values
    t850 = data01['t'].values
    z500 = data03['z'].values/9.81
    z850 = data01['z'].values/9.81
    lon = data01["longitude"]
    lat = data01["latitude"]
    X,Y = np.meshgrid(lon,lat)

    lapse_rate = (t500-t850)/(z500-z850)*(-1)
    
    # Set figure size to match geographic extent
    lat_extent = latN - latS
    lon_extent = lonE - lonW
    aspect_ratio = lon_extent / lat_extent

    base_height = 20  # inches
    fig_width = base_height * aspect_ratio  # 15 * 1 = 15

    fig,ax = plt.subplots(figsize=(fig_width,base_height), subplot_kw={'projection':ccrs.PlateCarree()},
                          constrained_layout=True)

#    fig, ax = plt.subplots(figsize=(20, 15),subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle='-')
    ax.add_feature(cfeature.STATES, linestyle='-')
    ax.add_feature(cfeature.COASTLINE)
    ax.gridlines(draw_labels=True)
    
    vmin = 0  #np.min(q)
    vmax = 20 #np.max(q)
    levels = 20
    level_boundaries = np.linspace(vmin, vmax, levels + 1)

    cf = ax.contourf(X, Y, q500*1000, level_boundaries, cmap='YlGn', alpha=1, vmin=0, vmax=20)


    shrinkscale=0.55
    
    cbar=fig.colorbar(
        ScalarMappable(norm=cf.norm, cmap=cf.cmap),
        ticks=range(vmin, vmax+5, 5),
        boundaries=level_boundaries,
        values=(level_boundaries[:-1] + level_boundaries[1:])/2,
        ax=ax, orientation='vertical', shrink=0.7, pad=0.05
    )
    cbar.set_label('Specific humidity [g/kg]', fontsize=16)  
    cbar.ax.tick_params(labelsize=12) 

    ct = ax.contour(X, Y, q500*1000, levels=[5,10,12,14,16,18,20], alpha=0.7, colors='white', linewidths=1)
    ax.clabel(ct, inline=True, fontsize=11, fmt='%d')

    ct8 = ax.contour(X, Y, lapse_rate*1000, levels=[8.0], colors='purple', linewidths=1)
    ax.clabel(ct8, inline=True, fontsize=11, fmt='%.1f')

    ct7 = ax.contour(X, Y, lapse_rate*1000, levels=[7.0], colors='red', linewidths=1)
    ax.clabel(ct7, inline=True, fontsize=11, fmt='%.1f')

    ct65 = ax.contour(X, Y, lapse_rate*1000, levels=[6.5], colors='red', linestyles='dashed', linewidths=1)
    ax.clabel(ct65, inline=True, fontsize=11, fmt='%.1f')

    ct60 = ax.contour(X, Y, lapse_rate*1000, levels=[6.], colors='red', linestyles='dashed', linewidths=0.5)
    ax.clabel(ct60, inline=True, fontsize=11, fmt='%.1f')
    
    #Wind arrows are plotted every n gridpoint
    cuv = ax.quiver(X[::n,::n], Y[::n,::n], u500[::n,::n], v500[::n,::n], angles='xy', scale_units='xy', scale=15)
    ax.quiverkey(cuv, X=1.09, Y=0.99, U=5, label='5 m/s', labelpos='E')
    
    plt.title('Specific humidity at 500m AGL (g/kg, shaded+white), Lapse Rate (K/km, red)\n\
    %sUTC' %(timt), fontsize = 20)
#    plt.text(lontstart,lattstart,"T", fontsize = 18, weight="bold",color='yellow')
    if case == 'tornado':
        ax.scatter(lontstart, lattstart, s=250, marker="v", edgecolors='k', color='yellow')
        ax.scatter(lont, latt, s=250, marker="v", edgecolors='k', color='green', alpha=0.7)
    else:
        ax.scatter(lontstart, lattstart, s=250, marker="^", edgecolors='k', color='cyan')
        ax.scatter(lont, latt, s=250, marker="^", edgecolors='k', color='green', alpha=0.7)
        
#    plt.tight_layout()
#    plt.savefig(pngfolder+'moisture_map.png')
    if pngfolder is not None:
        plt.savefig(pngfolder+"moisture_map.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
    plt.show()
    
    return u500,v500,q500,lapse_rate
    



def PlotCAPEMap(data,datas,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,pngfolder=None):
    data01 = data.sel(pressure_level=500, valid_time=timt)
#    data03 = data.sel(pressure_level=925, valid_time=timt)
#    data04 = data.sel(pressure_level=700, valid_time=timt)
#    data02 = datas.sel(valid_time=timt)
    data00 = data.sel(valid_time=timt)
    data00s = datas.sel(valid_time=timt)


#    q = data03["q"] 
    z = data01["z"]/9.81
#    uu = data03["u"]
#    vv = data03["v"]
#    uu2 = data04["u"]
#    vv2 = data04["v"]
    cape = data00s["cape"]
    u10s = data00s['u10'].values*1.94384449 #in knots
    v10s = data00s['v10'].values*1.94384449 #in knots
    orog = data00s['z'].values/9.81
    u00 = data00['u'].values*1.94384449 #in knots
    v00 = data00['v'].values*1.94384449 #in knots
    z00 = data00['z'].values/9.81
    q00 = data00['q'].values
    
    u500 = InterpolateToHeightAboveGround(z00, orog, np.array(u00), 500)
    v500 = InterpolateToHeightAboveGround(z00, orog, np.array(v00), 500)
    q500 = InterpolateToHeightAboveGround(z00, orog, np.array(q00), 500)

    u3000 = InterpolateToHeightAboveGround(z00, orog, np.array(u00), 3000)
    v3000 = InterpolateToHeightAboveGround(z00, orog, np.array(v00), 3000)
    
    lat = data['latitude']
    lon = data['longitude']
    X,Y = np.meshgrid(lon,lat)
    
    #    color_cols = colorscales['cape']['colors']
#    color_vals = colorscales['cape']['values']
#    mymap, mymin, mymax = Colormap(color_cols, color_vals)
    
    fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()})
    ax.set_extent([lonW,lonE,latS,latN], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle='-')
    ax.add_feature(cfeature.STATES, linestyle='-')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAKES, alpha=0.2)
    
    vmin = 100    #np.min(cape)
    vmax = 3000 #np.max(cape)
    nlevels = 32
    #colormap=cm.get_cmap("plasma").copy()
    # colormap = cm.get_cmap("inferno").copy()
    # colormap = pyart.graph.cmweather.cm.LangRainbow12.copy()
    colormap = cmocean.cm.thermal.copy()
    colormap.set_over("gray")
    colormap.set_under(alpha=0)

    alphas = np.maximum(0, np.minimum(1, (np.array(cape)/100)))
    alphas[np.isnan(alphas)] = 0
    
    # print(np.shape(X), np.shape(alphas))
    
    if Y[1,1] > Y[0,0]:
        myorigin = "lower"
    else:
        myorigin = "upper"

    if region == 'Canada':
        shrinkscale = 0.3
    else:
        shrinkscale = 0.6

    #cape2=np.array(cape)
    #cape2[cape2<100]=np.nan
    levels = np.arange(100,3001,100)
    
    alphas = np.zeros(len(levels))
    for i in range(0,len(levels)):
        alphas[i] = 1 - i/len(levels)
    norm = BoundaryNorm(np.linspace(vmin, vmax, nlevels+1), ncolors=colormap.N)
    cf = ax.contourf(X, Y, cape, levels=levels, cmap=colormap, norm=norm, vmin=vmin, vmax=vmax, extend='both')
    cbar = plt.colorbar(cf, ax=ax, orientation='vertical', shrink=shrinkscale)
    cbar.set_label('CAPE [J/kg]', fontsize=16)  
    cbar.ax.tick_params(labelsize=12)
    
    #extent=[np.min(X),np.max(X),np.min(Y),np.max(Y)]
    #cape_shaded = ax.imshow(cape, interpolation='bilinear', transform=ccrs.PlateCarree(), origin=myorigin, 
    #                        cmap='plasma', vmin=vmin, vmax=vmax, alpha=alphas, 
    ##                        cmap=mymap, vmin=mymin, vmax=mymax, alpha=alphas, 
    #                        extent=extent, aspect='auto', zorder=20)
    #cape_shaded_colorbar  = ax.imshow(cape, transform=ccrs.PlateCarree(), origin=myorigin, cmap='plasma', 
    #                                  vmin=vmin, vmax=vmax, extent=[0,0.5,0,0.5], aspect='auto', zorder=20)
#    fig.colorbar(cape_shaded, ax=ax, orientation='vertical', shrink=shrinkscale)
#    cbar = plt.colorbar(cape_shaded, ax=ax, orientation='vertical', shrink=shrinkscale)

    ct = ax.contour(X, Y, q500*1000, levels=[8,10,12,14,16,18,20], colors='cyan', linewidths=2)
    ax.clabel(ct, inline=True, fontsize=11, fmt='%d')

    #n=7
    #Wind barbs are plotted every n gridpoint
    cuv = ax.barbs(X[::n,::n], Y[::n,::n], u3000[::n,::n]-u10s[::n,::n], v3000[::n,::n]-v10s[::n,::n], color="k",)

    cs = ax.contour(X, Y, z, levels=np.arange(4800,6501,50), colors='k', linewidths=2)#'white', linewidths=1)
    ax.clabel(cs, inline=True, fontsize=11, fmt='%d')
    
    plt.title('10m-3km wind difference (barbs, knots), 500 hPa Geopotential (gpdm, black contours),\n\
               CAPE (shaded), specific humidity 500m AGL (g/kg, cyan contours)\n %sUTC' %(timt[:13]), fontsize = 20)
    if case == 'tornado':
        ax.scatter(lontstart, lattstart, s=200, marker="v", edgecolors='k', color='yellow', linewidth=1.25)
        ax.scatter(lont, latt, s=200, marker="v", edgecolors='k', color='green', linewidth=1.25)
    else:
        ax.scatter(lontstart, lattstart, s=200, marker="^", edgecolors='k', color='cyan', linewidth=1.25)
        ax.scatter(lont, latt, s=200, marker="^", edgecolors='k', color='green', linewidth=1.25)

    if pngfolder is not None:
        plt.savefig(pngfolder+"cape_map.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
    plt.show()
    
    return




def PlotPressureMaps(data,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,pngfolder=None):
    for presl in [850,700,500,300]:
        data01 = data.sel(pressure_level=presl, valid_time=timt)

        t = data01["t"].values-273.15
        z = data01["z"].values/9.81
        uu = data01["u"]*1.94384449 #in knots
        vv = data01["v"]*1.94384449 #in knots
        vo = data01["vo"].values
        pv = data01["pv"].values
        lon = data01["longitude"].values
        lat = data01["latitude"].values
        X,Y = np.meshgrid(lon,lat)

        fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()})
        ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.BORDERS, linestyle='-')
        ax.add_feature(cfeature.STATES, linestyle='-')
        ax.add_feature(cfeature.COASTLINE)
        
        if (presl == 850):
            # color_vals = [-40,-32,-28,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,\
            #               8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42]
            color_vals = np.arange(-20, 42, 2)
            color_cols = ['#ffffff','#ffe0ff','#ffc1ff','#ffa2ff','#ef7df4','#d854e6','#c22cd7','#a621d6',\
                          '#823ae4','#5f52f3','#416cfe','#3f8bf8','#3eaaf3','#3cc9ed','#51cdc0','#69cf8f',\
                          '#80d15d','#9dd947','#bee341','#deed3b','#f4ec36','#f8d533','#fcbe30','#ffa62c',\
                          '#ff8c1e','#ff7210','#ff5703','#ff3e36','#ff2579','#ff0bbb','#ff17df','#ff41eb',\
                          '#ff6af7','#ff92ff','#ffb7ff','#ffdbff','#ffffff']
            xmin_z = 800; xmax_z = 1800; xint_z = 25     # z interval
        
        elif (presl == 700):
            # color_vals = [-50,-40,-38,-36,-34,-32,-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,\
            #               4,6,8,10,12,14,16,18,20,22,30]
            color_vals = np.arange(-30, 32, 2)
            color_cols = ['#ffffff','#ffe0ff','#ffc1ff','#ffa2ff','#ef7df4','#d854e6','#c22cd7','#a621d6',\
                          '#823ae4','#5f52f3','#416cfe','#3f8bf8','#3eaaf3','#3cc9ed','#51cdc0','#69cf8f',\
                          '#80d15d','#9dd947','#bee341','#deed3b','#f4ec36','#f8d533','#fcbe30','#ffa62c',\
                          '#ff8c1e','#ff7210','#ff5703','#ff3e36','#ff2579','#ff0bbb','#ff17df','#ff41eb',\
                          '#ff6af7','#ff92ff','#ffb7ff','#ffdbff','#ffffff']
            xmin_z = 2400; xmax_z = 4000; xint_z = 50    # z interval

        elif (presl == 500):
            # color_vals = [-60,-56,-54,-52,-48,-46,-44,-42,-40,-38,-36,-34,-32,-30,-28,-26,-24,-22,\
            #               -20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,10]
            color_vals = np.arange(-40, 22, 2)
            color_cols = ['#ffffff','#ffe0ff','#ffc1ff','#ffa2ff','#ef7df4','#d854e6','#c22cd7','#a621d6',\
                          '#823ae4','#5f52f3','#416cfe','#3f8bf8','#3eaaf3','#3cc9ed','#51cdc0','#69cf8f',\
                          '#80d15d','#9dd947','#bee341','#deed3b','#f4ec36','#f8d533','#fcbe30','#ffa62c',\
                          '#ff8c1e','#ff7210','#ff5703','#ff3e36','#ff2579','#ff0bbb','#ff17df','#ff41eb',\
                          '#ff6af7','#ff92ff','#ffb7ff','#ffdbff','#ffffff']
            xmin_z = 4800; xmax_z = 6500; xint_z = 50    # z interval
        
        elif (presl == 300):
            # color_vals = [-90,-82,-80,-78,-76,-74,-72,-70,-68,-66,-64,-62,-60,-58,-56,-54,-52,-48,\
            #               -46,-44,-42,-40,-38,-36,-34,-32,-30,-28,-26,-24,-22,-20,-18,-16,-14,-10]
            color_vals = np.arange(-70, -8, 2)
            color_cols = ['#ffffff','#ffe0ff','#ffc1ff','#ffa2ff','#ef7df4','#d854e6','#c22cd7','#a621d6',\
                          '#823ae4','#5f52f3','#416cfe','#3f8bf8','#3eaaf3','#3cc9ed','#51cdc0','#69cf8f',\
                          '#80d15d','#9dd947','#bee341','#deed3b','#f4ec36','#f8d533','#fcbe30','#ffa62c',\
                          '#ff8c1e','#ff7210','#ff5703','#ff3e36','#ff2579','#ff0bbb','#ff17df','#ff41eb',\
                          '#ff6af7','#ff92ff','#ffb7ff','#ffdbff','#ffffff']
            xmin_z = 8000; xmax_z = 10000; xint_z = 100# z interval

        vmin = np.min(np.min(color_vals))#levels)#-20#np.min(t)
        vmax = np.max(np.max(color_vals)) #-5#np.max(t)

        if region == 'Canada':
            shrinkscale = 0.3
        else:
            shrinkscale = 0.6
        
        # cf = ax.contourf(X, Y, t, color_vals, colors=color_cols, alpha=1, vmin=vmin, vmax=vmax)
        cf = ax.contourf(X, Y, t, color_vals, cmap='LangRainbow12', alpha=1, vmin=vmin, vmax=vmax)
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Specific humidity [°C]', shrink=shrinkscale)
        cbar.set_label('Temperature [°C]', fontsize=16)  
        cbar.ax.tick_params(labelsize=12) 

        ct = ax.contour(X, Y, t, levels=30, colors='white', linewidths=1)
        ax.clabel(ct, inline=True, fontsize=11, fmt='%d')#'%.2f')
        if (presl == 850) | (presl == 700):
            c0 = ax.contour(X, Y, t, levels=[0], colors='b', linewidths=1.5, linestyles='--')
            ax.clabel(c0, inline=True, fontsize=11, fmt='%d')
    
        cs = ax.contour(X, Y, z, levels=np.arange(xmin_z,xmax_z+1,xint_z), colors='black', linewidths=2)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')
    
        #Wind barbs are plotted every n gridpoint
        cuv = ax.barbs(X[::n,::n], Y[::n,::n], uu[::n,::n], vv[::n,::n])
        
        if (presl == 300):
            windmag = np.sqrt((uu/1.94384449)**2 + (vv/1.94384449)**2)
            windmag = np.array(windmag)
            windmag[windmag<30] = np.nan
            windmag[windmag>=30] = 1
            ax.contourf(X, Y, windmag, levels=[1,1.1], colors='none', alpha=0.2)
            plt.title('%d hPa Geopotential [gpdm, contour], Temperature (color), Wind>30m/s (dark shaded)\n\
            %sUTC' %(presl,timt), fontsize=20)
        else:
            plt.title('%d hPa Geopotential [gpdm, contour], Temperature (color), Wind (barbs, knots)\n\
            %sUTC' %(presl,timt), fontsize=20)
        if case == 'tornado':
            ax.scatter(lontstart, lattstart, s=200, marker="v", edgecolors='k', color='magenta', linewidth=1.25)
            ax.scatter(lont, latt, s=200, marker="v", edgecolors='k', color='green', linewidth=1.25)
        else:
            ax.scatter(lontstart, lattstart, s=200, marker="^", edgecolors='k', color='cyan', linewidth=1.25)
            ax.scatter(lont, latt, s=200, marker="^", edgecolors='k', color='green', linewidth=1.25)
        
        # iyl,ixl = mc.find_peaks(z, maxima=False)
        # iyh,ixh = mc.find_peaks(z)
        # for i in range(len(iyl)):
        #     lonmin = lon[ixl[i]]
        #     latmin = lat[iyl[i]]
        #     ax.scatter(lonmin, latmin, s=400, marker='*', color='k', facecolor='b', linewidth=1.5)
        # for j in range(len(iyh)):
        #     lonmax = lon[ixh[j]]
        #     latmax = lat[iyh[j]]
        #     ax.scatter(lonmax, latmax, s=400, marker='*', color='k', facecolor='r', linewidth=1.5)
        
        
        if pngfolder is not None:
            plt.savefig(pngfolder+f"pressure_map_{presl}mb.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
        # plt.show()
        
    return




def PlotMoistureMap2(data,datas,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,pngfolder=None):
#    data01 = data.sel(pressure_level=925, valid_time=timt)
    data01 = data.sel(pressure_level=850, valid_time=timt)
    data03 = data.sel(pressure_level=500, valid_time=timt)
    data02 = datas.sel(valid_time=timt)
    data00 = data.sel(valid_time=timt)
    orog = data02['z'].values/9.81
    z00 = data00['z'].values/9.81
    u00 = data00['u'].values
    v00 = data00['v'].values
    z00 = data00['z'].values/9.81
    q00 = data00['q'].values
    
    u500 = InterpolateToHeightAboveGround(z00, orog, np.array(u00), 500)
    v500 = InterpolateToHeightAboveGround(z00, orog, np.array(v00), 500)
    q500 = InterpolateToHeightAboveGround(z00, orog, np.array(q00), 500)

#    q = data01["q"] 
#    z = data01["z"] / (9.81*10) #in gpdm
#    uu = data01["u"]
#    vv = data01["v"]
    t500 = data03['t'].values
    t850 = data01['t'].values
    z500 = data03['z'].values/9.81
    z850 = data01['z'].values/9.81
    lon = data01["longitude"]
    lat = data01["latitude"]
    X,Y = np.meshgrid(lon,lat)

    lapse_rate = (t500-t850)/(z500-z850)*(-1)
    
    fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()})
    ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle='-')
    ax.add_feature(cfeature.STATES, linestyle='-')
    ax.add_feature(cfeature.COASTLINE)
    
    vmin = 0  #np.min(q)
    vmax = 20 #np.max(q)
    levels = 20
    level_boundaries = np.linspace(vmin, vmax, levels+1)

    cf = ax.contourf(X, Y, q500*1000, level_boundaries, cmap='YlGn', alpha=1, vmin=0, vmax=20)


    if region == 'Canada':
        shrinkscale = 0.25
        scale = 5
    else:
        shrinkscale = 0.55
        scale = 10
        
    cbar = fig.colorbar(
        ScalarMappable(norm=cf.norm, cmap=cf.cmap),
        ticks=range(vmin, vmax+5, 5),
        boundaries=level_boundaries,
        values=(level_boundaries[:-1] + level_boundaries[1:]) / 2,
        ax=ax, orientation='vertical', shrink=shrinkscale
    )
    cbar.set_label('Specific humidity [g/kg]', fontsize=16)  
    cbar.ax.tick_params(labelsize=12) 

    ct = ax.contour(X, Y, q500*1000, levels=[5,10,12,14,16,18,20], alpha=0.7, colors='white', linewidths=1)
    ax.clabel(ct, inline=True, fontsize=11, fmt='%d')

    ct8 = ax.contour(X, Y, lapse_rate*1000, levels=[8.0], colors='purple', linewidths=1)
    ax.clabel(ct8, inline=True, fontsize=11, fmt='%.1f')

    ct7 = ax.contour(X, Y, lapse_rate*1000, levels=[7.0], colors='red', linewidths=1)
    ax.clabel(ct7, inline=True, fontsize=11, fmt='%.1f')

    ct65 = ax.contour(X, Y, lapse_rate*1000, levels=[6.5], colors='red', linestyles='dashed', linewidths=1)
    ax.clabel(ct65, inline=True, fontsize=11, fmt='%.1f')

    ct60 = ax.contour(X, Y, lapse_rate*1000, levels=[6.], colors='red', linestyles='dashed', linewidths=0.5)
    ax.clabel(ct60, inline=True, fontsize=11, fmt='%.1f')

    #Wind arrows are plotted every n gridpoint
    cuv = ax.quiver(X[::n,::n], Y[::n,::n], u500[::n,::n], v500[::n,::n], angles='xy', scale_units='xy', scale=scale)
    ax.quiverkey(cuv, X=1.08, Y=0.99, U=scale, label=f"{scale} m/s", labelpos='N', fontproperties=dict(size=12))
    
    plt.title('Specific humidity at 500m AGL (g/kg, shaded+white), Lapse Rate (K/km, red)\n\
    %sUTC' %(timt), fontsize=20)
#    plt.text(lontstart,lattstart,"T", fontsize = 18, weight="bold",color='yellow')
    if case == 'tornado':
        ax.scatter(lontstart, lattstart, s=200, marker="v", edgecolors='k', color='magenta', linewidth=1.25)
        ax.scatter(lont, latt, s=200, marker="v", edgecolors='k', color='green', linewidth=1.25)
    else:
        ax.scatter(lontstart, lattstart, s=200, marker="^", edgecolors='k', color='cyan', linewidth=1.25)
        ax.scatter(lont, latt, s=200, marker="^", edgecolors='k', color='green', linewidth=1.25)
    
    if pngfolder is not None:
        plt.savefig(pngfolder+"moisture_map.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
    plt.show()
    
    return




def PlotSRHMaps(data,datas,lonW,lonE,latS,latN,timt,lont,latt,lontstart,lattstart,region,case,n,pngfolder=None):
    
    data00 = data.sel(valid_time=timt)
    data00s = datas.sel(valid_time=timt)
    
    u10s = data00s['u10'].values
    v10s = data00s['v10'].values
    orog = data00s['z'].values/9.81
    u00 = data00['u']
    v00 = data00['v']
    z00 = data00['z'].values/9.81
    p00 = data00.pressure_level.values
    z500 = z00[(p00==500),:,:].squeeze()
    
    lat = data00['latitude']
    lon = data00['longitude']
    X,Y = np.meshgrid(lon,lat)
    
    
    
    for zdep in [1000,3000]:
        
        uu = InterpolateToHeightAboveGround(z00, orog, np.array(u00.values), zdep)
        vv = InterpolateToHeightAboveGround(z00, orog, np.array(v00.values), zdep)
        uu = uu - u10s
        vv = vv - v10s
        
        
        [srhpos,srhneg,srhtot] = getSRH(z00, orog, u00.values, v00.values, zdep, u10s, v10s, storm_vector='rightmover')
        
        srh = srhtot.magnitude
        
                
        
        fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()})
        ax.set_extent([lonW,lonE,latS,latN], crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.BORDERS, linestyle='-')
        ax.add_feature(cfeature.STATES, linestyle='-')
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.LAKES, alpha=0.2)
        
        if zdep == 3000:
            vmin = 10    #np.min(srh)
            vmax = 300 #np.max(srh)
            nlevels = 30
        elif zdep == 1000:
            vmin = 10    #np.min(srh)
            vmax = 300 #np.max(srh)
            nlevels = 30
        
        # vmin = 10
        # vmax = 1000
        # nlevels = 50
        
        
        colormap = pyart.graph.cmweather.cm.LangRainbow12.copy()
        colormap.set_over("gray")
        colormap.set_under(alpha=0)
    
        alphas = np.maximum(0, np.minimum(1, (np.array(srh)/10)))
        alphas[np.isnan(alphas)] = 0
        
        # print(np.shape(X), np.shape(alphas))
        
        if Y[1,1] > Y[0,0]:
            myorigin = "lower"
        else:
            myorigin = "upper"
    
        if region == 'Canada':
            shrinkscale = 0.3
            scale = 5
        else:
            shrinkscale = 0.6
            scale = 10
    
        
        if zdep == 3000:
            levels = np.arange(10,301,10)
        elif zdep == 1000:
            levels = np.arange(10,301,10)
        
        
        
        alphas = np.zeros(len(levels))
        for i in range(0,len(levels)):
            alphas[i] = 1 - i/len(levels)
        norm = BoundaryNorm(np.linspace(vmin, vmax, nlevels+1), ncolors=colormap.N)
        
        cf = ax.contourf(X, Y, srh, levels=levels, cmap=colormap, norm=norm, vmin=vmin, vmax=vmax, extend='both')
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', shrink=shrinkscale)
        cbar.set_label('SRH [m2/s2]', fontsize=16)  
        cbar.ax.tick_params(labelsize=12)
        
        
    
    
        #n=7
        #Wind barbs are plotted every n gridpoint
        cuv = ax.quiver(X[::n,::n], Y[::n,::n], uu[::n,::n], vv[::n,::n], angles='xy', scale_units='xy', scale=scale)
        ax.quiverkey(cuv, X=1.08, Y=1.02, U=10, label="10 m/s", labelpos='N', fontproperties=dict(size=12))
    
        cs = ax.contour(X, Y, z500, levels=np.arange(4800,6501,50), colors='k', linewidths=2)#'white', linewidths=1)
        ax.clabel(cs, inline=True, fontsize=11, fmt='%d')
        
        plt.title(f" 10m-{zdep/1000:.0f}km wind difference (vectors, m/s), 500 hPa Geopotential (gpdm, black contours),\n\
                  0-{zdep/1000:.0f}km SRH (shaded)\n {timt[:13]} UTC", fontsize=20)
        if case == 'tornado':
            ax.scatter(lontstart, lattstart, s=200, marker="v", edgecolors='k', color='yellow', linewidth=1.25)
            ax.scatter(lont, latt, s=200, marker="v", edgecolors='k', color='green', linewidth=1.25)
        else:
            ax.scatter(lontstart, lattstart, s=200, marker="^", edgecolors='k', color='cyan', linewidth=1.25)
            ax.scatter(lont, latt, s=200, marker="^", edgecolors='k', color='green', linewidth=1.25)
    
        if pngfolder is not None:
            plt.savefig(pngfolder+"srh_0{zdep/1000:.0f}km_map.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
        # plt.show()
    
    return




# Convert lat/lon coordinates to x/y distances relative to an origin point (in km)
def latlon2xy(lat, lon, lat_o, lon_o):
    # lat, lon:     1-D vectors of lat/lon in decimal degrees N/deg E
    # lat_o, lon_o: lat/lon of origin in decimal degrees N/deg E
    
    r_earth = 6378.1 # km
    
    thy = lat_o*np.pi/180 # convert to radians
    thz = -lon_o*np.pi/180
    
    # transform matrices
    Ry = [[np.cos(thy), 0, np.sin(thy)],
          [0, 1, 0],
          [-np.sin(thy), 0, np.cos(thy)]]
    
    Rz = [[np.cos(thz), -np.sin(thz), 0],
          [np.sin(thz), np.cos(thz), 0],
          [0, 0, 1]]
    
    # i'm actually not sure exactly how this works, i just copied this function from some of Boonleng's code
    R = np.matmul(Ry,Rz)
    xyz = r_earth * np.array([np.cos(lat*np.pi/180) * np.cos(lon*np.pi/180),
                np.cos(lat*np.pi/180) * np.sin(lon*np.pi/180),
                np.sin(lat*np.pi/180)])
    # get x and y positions
    posx = np.matmul(R[1],xyz)
    posy = np.matmul(R[2],xyz)
    
    # if len(lat) != len(lon):
    #     posx = np.zeros(shape=(len(lat),len(lon)))
    #     posy = np.zeros(shape=(len(lat),len(lon)))
    #     for i in range(len(lon)):
    #         for j in range(len(lat)):
    #             xyz = r_earth * np.array([np.cos(lat[j]*np.pi/180) * np.cos(lon[i]*np.pi/180),
    #                         np.cos(lat[j]*np.pi/180) * np.sin(lon[i]*np.pi/180),
    #                         np.sin(lat[j]*np.pi/180)])
    #             posx[j,i] = np.matmul(R[1],xyz)
    #             posy[j,i] = np.matmul(R[2],xyz)
    # else:
    #     xyz = r_earth * np.array([np.cos(lat*np.pi/180) * np.cos(lon*np.pi/180),
    #                 np.cos(lat*np.pi/180) * np.sin(lon*np.pi/180),
    #                 np.sin(lat*np.pi/180)])
    #     # get x and y positions
    #     posx = np.matmul(R[1],xyz)
    #     posy = np.matmul(R[2],xyz)
    
    return posx,posy






def getStormMotion(z, orog, u, v):
    # orog = data00s["z"].values/9.81
    # u = data00['u']
    # v = data00['v']
    # z = data00['z'].values/9.81
    
    u500 = InterpolateToHeightAboveGround(z, orog, u, 500)
    v500 = InterpolateToHeightAboveGround(z, orog, v, 500)
    u5500 = InterpolateToHeightAboveGround(z, orog, u, 5500)
    v5500 = InterpolateToHeightAboveGround(z, orog, v, 5500)
    u6000 = InterpolateToHeightAboveGround(z, orog, u, 6000)
    v6000 = InterpolateToHeightAboveGround(z, orog, v, 6000)
    
    
    z_interp = z
    u_interp = u
    v_interp = v
    
    if 6000 not in z_interp:
        arr1 = 6000 * np.ones(shape=orog.shape)
        z1 = np.append(z_interp, arr1[np.newaxis,:], axis=0)
        sort_inds1 = np.argsort(z1[:,0,0])
        z_interp = z1[sort_inds1,:,:]
        
        u1 = np.append(u_interp, u6000[np.newaxis,:], axis=0)
        u_interp = u1[sort_inds1,:,:]
        v1 = np.append(v_interp, v6000[np.newaxis,:], axis=0)
        v_interp = v1[sort_inds1,:,:]
        
    if 5500 not in z_interp:
        arr2 = 5500 * np.ones(shape=orog.shape)
        z2 = np.append(z_interp, arr2[np.newaxis,:], axis=0)
        sort_inds2 = np.argsort(z2[:,0,0])
        z_interp = z2[sort_inds2,:,:]
        
        u2 = np.append(u_interp, u5500[np.newaxis,:], axis=0)
        u_interp = u2[sort_inds2,:,:]
        v2 = np.append(v_interp, v5500[np.newaxis,:], axis=0)
        v_interp = v2[sort_inds2,:,:]
        
    if 500 not in z_interp:
        arr3 = 500 * np.ones(shape=orog.shape)
        z3 = np.append(z_interp, arr3[np.newaxis,:], axis=0)
        sort_inds3 = np.argsort(z3[:,0,0])
        z_interp = z3[sort_inds3,:,:]
        
        u3 = np.append(u_interp, u500[np.newaxis,:], axis=0)
        u_interp = u3[sort_inds3,:,:]
        v3 = np.append(v_interp, v500[np.newaxis,:], axis=0)
        v_interp = v3[sort_inds3,:,:]
    
    mask1 = (z_interp<6000) | np.isclose(z_interp,6000)
    u_0_6 = np.ma.masked_array(u_interp, ~mask1)
    v_0_6 = np.ma.masked_array(v_interp, ~mask1)
    
    mask2 = (z_interp<500) | np.isclose(z_interp,500)
    u_0_05 = np.ma.masked_array(u_interp, ~mask2)
    v_0_05 = np.ma.masked_array(v_interp, ~mask2)
    
    mask3 = ((z_interp>5500) | np.isclose(z_interp,5500)) & ((z_interp<6000) | np.isclose(z_interp,6000))
    u_55_6 = np.ma.masked_array(u_interp, ~mask3)
    v_55_6 = np.ma.masked_array(v_interp, ~mask3)
    
    u_mean = np.nanmean(u_0_6, axis=0)
    v_mean = np.nanmean(v_0_6, axis=0)
    wind_mean = np.array([u_mean, v_mean])
    
    u_500m = np.nanmean(u_0_05, axis=0)
    v_500m = np.nanmean(v_0_05, axis=0)
    wind_500m = np.array([u_500m, v_500m])
    
    u_5500m = np.nanmean(u_55_6, axis=0)
    v_5500m = np.nanmean(v_55_6, axis=0)
    wind_5500m = np.array([u_5500m, v_5500m])
    
    shear = wind_5500m - wind_500m
    shear_cross = np.asarray([shear[1,:,:], -1*shear[0,:,:]])
    shear_mag = np.sqrt(shear[0,:,:]**2 + shear[1,:,:]**2)
    rdev = shear_cross * (7.5 / shear_mag)
    
    right_mover = wind_mean + rdev
    left_mover = wind_mean - rdev
    
    return right_mover, left_mover, wind_mean




def getSRH(z, orog, u, v, depth, u10, v10, bottom=None, storm_vector=None):
    
    if bottom is None:
        bottom = 10
        top = depth
    # if storm_u is None:
    #     storm_u = 0
    # if storm_v is None:
    #     storm_v = 0
    else:
        top = bottom + depth
    
    z_interp = z
    u_interp = u
    v_interp = v
    
    if top not in z_interp:
        top_arr = top * np.ones(shape=orog.shape)
        z1 = np.append(z_interp, top_arr[np.newaxis,:], axis=0)
        sort_inds1 = np.argsort(z1[:,0,0])
        
        z_interp = z1[sort_inds1,:,:]
        
        u_top = InterpolateToHeightAboveGround(z, orog, u, top)
        u1 = np.append(u_interp, u_top[np.newaxis,:], axis=0)
        u_interp = u1[sort_inds1,:,:]

        v_top = InterpolateToHeightAboveGround(z, orog, v, top)
        v1 = np.append(v_interp, v_top[np.newaxis,:], axis=0)
        v_interp = v1[sort_inds1,:,:]
        
    # if bottom not in z_interp:
    if True:
        bot_arr = bottom * np.ones(shape=orog.shape)
        z2 = np.append(z_interp, bot_arr[np.newaxis,:], axis=0)
        sort_inds2 = np.argsort(z2[:,0,0])
        
        z_interp = z2[sort_inds2,:,:]
        
        u2 = np.append(u_interp, u10[np.newaxis,:], axis=0)
        u_interp = u2[sort_inds2,:,:]
        
        v2 = np.append(v_interp, v10[np.newaxis,:], axis=0)
        v_interp = v2[sort_inds2,:,:]
    
    
    [rm,lm,mw] = getStormMotion(z, orog, u, v)
    u_rm=rm[0,:,:]; v_rm=rm[1,:,:]
    u_lm=lm[0,:,:]; v_lm=lm[1,:,:]
    u_mw=mw[0,:,:]; v_mw=mw[1,:,:]
    
    if storm_vector is None:
        storm_u = u_rm
        storm_v = v_rm
    elif (storm_vector == 'rightmover') | (storm_vector == 'right'):
        storm_u = u_rm
        storm_v = v_rm
    elif (storm_vector == 'leftmover') | (storm_vector == 'left'):
        storm_u = u_lm
        storm_v = v_lm
    elif (storm_vector == 'meanwind') | (storm_vector == 'mean'):
        storm_u = u_mw
        storm_v = v_mw
    
    mask = ((z_interp>bottom) | np.isclose(z_interp,bottom)) & ((z_interp<top) | np.isclose(z_interp,top))
    u_layer = np.ma.masked_array(u_interp, ~mask)
    v_layer = np.ma.masked_array(v_interp, ~mask)
    
    storm_relative_u = u_layer - np.tile(storm_u, [u_layer.shape[0], 1, 1])
    storm_relative_v = v_layer - np.tile(storm_v, [v_layer.shape[0], 1, 1])
    
    int_layers = (storm_relative_u[1:,:,:] * storm_relative_v[:-1,:,:]
                  - storm_relative_u[:-1,:,:] * storm_relative_v[1:,:,:])
    int_layers_pos = np.ma.masked_array(int_layers, int_layers<0)
    int_layers_neg = np.ma.masked_array(int_layers, int_layers>0)
    
    positive_srh = np.asarray(np.nansum(int_layers_pos, axis=0))
    negative_srh = np.asarray(np.nansum(int_layers_neg, axis=0))
    total_srh = np.asarray(positive_srh + negative_srh)
    
    return (positive_srh*units('meter**2 / second**2'),
            negative_srh*units('meter**2 / second**2'),
            total_srh*units('meter**2 / second**2'))




















