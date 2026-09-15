# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 11:15:55 2026

@author: mschne28
"""

# Calculate QTor

import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import xarray as xr
import pandas as pd
import pyart #need an earlier version of xarray -> 0.20.2 or earlier
import pickle
import metpy.calc as mc
from metpy.plots import SkewT, Hodograph
from metpy.units import units
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from glob import glob
from scipy.ndimage import gaussian_filter
import os
from os.path import exists
from matplotlib.ticker import MultipleLocator


### STEPS
# 1. QLCS object identification DONEEEEEE
# 2. QLCS leading line identification
# 3. Storm motion estimation
# 4. Leading line shear analysis
# 5. Local tortuosity
# 6. Calculate QTor


# NEXRAD fields: 'differential_reflectivity', 'reflectivity', 'cross_correlation_ratio', 'clutter_filter_power_removed', 'spectrum_width', 'differential_phase', 'velocity'



#%% QLCS object criteria
# Britt et al. 2026 https://doi.org/10.1175/WAF-D-25-0092.1

# Britt et al. 2024 QLCS criteria https://doi.org/10.1175/WAF-D-23-0106.1
min_cref = 40 #composite dbz threshold to get initial objects
max_cref = 45 #max object dbz threshold
min_area = 54 #storm object area threshold [km^2]
merge_dist = 12 #proximity threshold [km]
min_length_1 = 100 #length threshold [km]
min_length_2 = 150 #length threshold [km]
min_ecc_1 = 0.85 #eccentricity threshold for length>100 km
min_ecc_2 = 0.74 #eccentricity threshold for length>150 km








#%% Open files

fp = 'C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/qlcs_tornado_outbreaks/'

radar_name = 'KBUF'
yyyyt = 2026
mmt = 8
ddt = 2
timestr = '162824'
filename = fp + f"{radar_name}/{radar_name}{yyyyt}{mmt:02.0f}{ddt:02.0f}_{timestr}_V06.ar2v"


radar = pyart.io.read(filename)

gate_x = radar.extract_sweeps([0]).gate_x['data']/1000
gate_y = radar.extract_sweeps([0]).gate_y['data']/1000
rng = radar.range['data']/1000

# dbz = radar.fields['reflectivity']['data']
gatefilter = pyart.filters.GateFilter(radar)
gatefilter.exclude_transition()
gatefilter.exclude_below('cross_correlation_ratio', 0.9)
cref = pyart.retrieve.composite_reflectivity(radar, field='reflectivity', gatefilter=gatefilter).fields['composite_reflectivity']['data']

#%% Identify objects

from skimage import measure,morphology,filters
from CM1utils import *
from metpy.interpolate import interpolate_to_points


# turn this into a function



pts = np.array([gate_x[:,(rng<=300)].ravel(), gate_y[:,(rng<=300)].ravel()]).transpose()
vals = cref[:,(rng<=300)].ravel().data
xm,ym = np.meshgrid(np.arange(-300,303,3), np.arange(-300,303,3))
xi = np.array([xm.ravel(), ym.ravel()]).transpose()

cref_interp_flat = interpolate_to_points(pts, vals, xi, interp_type='linear', search_radius=3)
cref_interp = cref_interp_flat.reshape(xm.shape[1], xm.shape[0])
#%%

cref_smooth = filters.gaussian(cref_interp, sigma=0.8)



# cref_bin = np.zeros(shape=cref.shape, dtype=bool)
# cref_bin[(cref>41)] = True
# cref_bin_clean = morphology.remove_small_objects(cref_bin, min_size=864, connectivity=1)
cref_bin = np.zeros(shape=cref_smooth.shape, dtype=bool)
cref_bin[(cref_smooth>40)] = True
cref_bin_clean = morphology.remove_small_objects(cref_bin, min_size=6, connectivity=1)


cref_labeled = measure.label(cref_bin_clean, connectivity=2)

regions = measure.regionprops(cref_labeled)


gate_x = radar.extract_sweeps([0]).gate_x['data']/1000
gate_y = radar.extract_sweeps([0]).gate_y['data']/1000

# bbox = regions[3].bbox
# ix1 = bbox[0]
# iy1 = bbox[1]
# ix2 = bbox[2]-1
# iy2 = bbox[3]-1


fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(gate_x, gate_y, cref, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
ax.set_title('Original composite reflectivity')
plt.show()


fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_interp, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
ax.set_title('Composite reflectivity interpolated to 3-km grid')
plt.show()


fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
ax.set_title('Gridded CRef with Gaussian smoothing')
plt.show()

gate_x = xm
gate_y = ym


fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(gate_x, gate_y, cref_bin, 'dbz', ax, datalims=[0,1], cmap='Grays')
ax.set_title('Binary CRef > 40 dBZ')
plt.show()


fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(gate_x, gate_y, cref_bin_clean, 'dbz', ax, datalims=[0,1], cmap='Grays')
ax.set_title('Cleaned binary CRef > 40 dBZ')
plt.show()


fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(gate_x, gate_y, cref_labeled, 'dbz', ax, datalims=[0,len(regions)], cmap='HomeyerRainbow')
# ax.scatter([gate_x[ix1,iy1], gate_x[ix1,iy2], gate_x[ix2,iy1], gate_x[ix2,iy2]], [gate_y[ix1,iy1], gate_y[ix1,iy2], gate_y[ix2,iy1], gate_y[ix2,iy2]], marker='.', s=30, c='k')
ax.set_title('CRef objects')
plt.show()


#%%

regions_filtered = []
cref_labeled_filtered = np.zeros(shape=cref_labeled.shape, dtype=int)

for i in range(len(regions)):
    area = regions[i].area
    axis_major = regions[i].axis_major_length
    coords = regions[i].coords
    bbox = regions[i].bbox
    ecc = regions[i].eccentricity
    fdiam = regions[i].feret_diameter_max
    
    # print(f"Region {i} max cref: {np.nanmax(cref_smooth[(cref_labeled==i)]):.1f}")
    print(f"Region {i+1} major axis, ecc: {axis_major:.1f}, {ecc:.2f}")
    
    
    # maxz_condition = np.max(cref_smooth[(cref_labeled==i)]) < 45
    # area_condition = area < 15
    # length_condition = axis_major < 33
    
    if np.nanmax(cref_smooth[(cref_labeled==i+1)]) > 45:
        maxz_met = True
        cref_labeled_filtered[(cref_labeled==i+1)] = i+1
        regions_filtered.append(i+1)
    else:
        maxz_met = False
    
    # if (axis_major > 33) & (ecc > 0.85):
    #     lenecc_met = True
    # elif (axis_major > 50) & (ecc > 0.74):
    #     lenecc_met = True
    # else:
    #     lenecc_met = False
    
    # if maxz_met & lenecc_met:
    #     isQLCS = True
    # else:
    #     isQLCS = False
    
    
    # if isQLCS:
    #     cref_labeled_filtered[(cref_labeled==i)] = i
    
    
# regions_filtered = regions_filtered[(regions_filtered is not None)]




fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(gate_x, gate_y, cref_labeled_filtered, 'dbz', ax, datalims=[0,len(regions)], cmap='HomeyerRainbow')
ax.set_title('Filtered CRef objects with max CRef > 45 dBZ')
plt.show()

#%%

cref_bin_filtered = np.zeros(shape=cref_smooth.shape, dtype=bool)
cref_bin_filtered[(cref_labeled_filtered>0)] = True

footprint = morphology.footprint_rectangle((7,7))
cref_bin_dilated = morphology.binary_dilation(cref_bin_filtered, footprint=footprint)

cref_merged = measure.label(cref_bin_dilated, connectivity=2)
regions_merged = measure.regionprops(cref_merged)





fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(gate_x, gate_y, cref_bin_dilated, 'dbz', ax, datalims=[0,1], cmap='Grays')
ax.set_title('Dilated binary filtered CRef objects')
plt.show()


fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(gate_x, gate_y, cref_merged, 'dbz', ax, datalims=[0,len(regions_filtered)], cmap='HomeyerRainbow')
ax.set_title('Merged filtered CRef objects')
plt.show()



regions_fin = []
cref_merged_filtered = np.zeros(shape=cref_merged.shape, dtype=int)

for i in range(len(regions_merged)):
    area = regions_merged[i].area
    axis_major = regions_merged[i].axis_major_length
    ecc = regions_merged[i].eccentricity
    
    print(f"Region {i+1} major axis, ecc: {axis_major:.1f}, {ecc:.2f}")
    
    
    if (axis_major > 33) & (ecc > 0.85):
        lenecc_met = True
    elif (axis_major > 50) & (ecc > 0.74):
        lenecc_met = True
    else:
        lenecc_met = False
    
    if lenecc_met:
        isQLCS = True
    else:
        isQLCS = False
    
    
    if isQLCS:
        cref_merged_filtered[(cref_merged==i+1)] = i+1
        cref_merged_filtered[(cref_bin_filtered==0)] = 0
        regions_fin.append(i+1)
        


### YESSSS IT WORKS

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(gate_x, gate_y, cref_merged_filtered, 'dbz', ax, datalims=[0,len(regions_fin)], cmap='HomeyerRainbow')
ax.set_title('Final QLCS objects\n filtered for max dBZ, merged, and filtered for length and eccentricity ')
plt.show()







