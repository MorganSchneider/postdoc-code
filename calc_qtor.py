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

from skimage import measure,morphology,filters
from QTORutils import *
from metpy.interpolate import interpolate_to_points



### STEPS
# 1. QLCS object identification DONEEEEEE
# 2. QLCS leading line identification
# 3. Storm motion estimation
# 4. Leading line shear analysis
# 5. Local tortuosity
# 6. Calculate QTor


# NEXRAD fields: 'differential_reflectivity', 'reflectivity', 'cross_correlation_ratio', 'clutter_filter_power_removed', 'spectrum_width', 'differential_phase', 'velocity'











#%% Open files

fp = 'C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/qlcs_tornado_outbreaks/'

radar_name = 'KBUF'
yyyyt = 2026
mmt = 8
ddt = 2
# timestr = '162824'
timestr = '161342'
filename = fp + f"{radar_name}/{radar_name}{yyyyt}{mmt:02.0f}{ddt:02.0f}_{timestr}_V06.ar2v"


radar = pyart.io.read(filename)
# dbz = radar.fields['reflectivity']['data']
gatefilter = pyart.filters.GateFilter(radar)
gatefilter.exclude_transition()
gatefilter.exclude_below('cross_correlation_ratio', 0.9)
cradar = pyart.retrieve.composite_reflectivity(radar, field='reflectivity', gatefilter=gatefilter)
cref = cradar.fields['composite_reflectivity']['data']

gate_x = cradar.extract_sweeps([0]).gate_x['data']/1000
gate_y = cradar.extract_sweeps([0]).gate_y['data']/1000
rng = cradar.range['data']/1000

max_rng = 300 #set max range, km
hres = 3 #gridded resolution, km

pts = np.array([gate_x[:,(rng<=max_rng)].ravel(), gate_y[:,(rng<=max_rng)].ravel()]).transpose()
vals = cref[:,(rng<=max_rng)].ravel().data
xm,ym = np.meshgrid(np.arange(-max_rng, max_rng+hres, hres), np.arange(-max_rng, max_rng+hres, hres))
xi = np.array([xm.ravel(), ym.ravel()]).transpose()

cref_interp_flat = interpolate_to_points(pts, vals, xi, interp_type='linear', search_radius=hres)
cref_interp = cref_interp_flat.reshape(xm.shape)

cref_smooth = filters.gaussian(cref_interp, sigma=0.8)


#%% Retrieve QLCS object IDs


qlcs_labels, qlcs_regions = find_qlcs_objects(cref_smooth, hres, merge_distance=6)

# if more than one object, pick the one with the highest dbz
if len(qlcs_regions) > 1:
    crefmax = np.array([np.nanmax(cref_smooth[(qlcs_labels==i+1)]) for i in range(len(qlcs_regions))])
    qlcs_obj = np.zeros(shape=qlcs_labels.shape, dtype=int)
    qlcs_obj[(qlcs_labels == np.argmax(crefmax)+1)] = 1
    qlcs_region = qlcs_regions[np.argmax(crefmax)]
else:
    qlcs_obj = qlcs_labels
    qlcs_region = qlcs_regions[0]





# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, qlcs_labels, 'dbz', ax, datalims=[0,len(qlcs_regions)], cmap='HomeyerRainbow')
# ax.set_title('Final QLCS objects\n filtered for max dBZ, merged, and filtered for length and eccentricity ')
# plt.show()



#%% Storm motion estimation - get data from previous and next volumes


filenames = glob(fp+f"{radar_name}/{radar_name}{yyyyt}{mmt:02.0f}{ddt:02.0f}_*_V06.ar2v")
ind = [i for i in range(len(filenames)) if timestr in filenames[i]][0]


# Previous volume
radar1 = pyart.io.read(filenames[ind-2])
gatefilter = pyart.filters.GateFilter(radar1)
gatefilter.exclude_transition()
gatefilter.exclude_below('cross_correlation_ratio', 0.9)
cradar_prev = pyart.retrieve.composite_reflectivity(radar1, field='reflectivity', gatefilter=gatefilter)
cref_prev = cradar_prev.fields['composite_reflectivity']['data']

gate_x = cradar_prev.extract_sweeps([0]).gate_x['data']/1000
gate_y = cradar_prev.extract_sweeps([0]).gate_y['data']/1000

pts = np.array([gate_x[:,(rng<=max_rng)].ravel(), gate_y[:,(rng<=max_rng)].ravel()]).transpose()
vals = cref_prev[:,(rng<=max_rng)].ravel().data

cref_interp_flat = interpolate_to_points(pts, vals, xi, interp_type='linear', search_radius=hres)
cref_interp_prev = cref_interp_flat.reshape(xm.shape)
cref_smooth_prev = filters.gaussian(cref_interp_prev, sigma=0.8)




# Next volume
radar2 = pyart.io.read(filenames[ind+2])
gatefilter = pyart.filters.GateFilter(radar2)
gatefilter.exclude_transition()
gatefilter.exclude_below('cross_correlation_ratio', 0.9)
cradar_next = pyart.retrieve.composite_reflectivity(radar2, field='reflectivity', gatefilter=gatefilter)
cref_next = cradar_next.fields['composite_reflectivity']['data']

gate_x = cradar_next.extract_sweeps([0]).gate_x['data']/1000
gate_y = cradar_next.extract_sweeps([0]).gate_y['data']/1000

pts = np.array([gate_x[:,(rng<=max_rng)].ravel(), gate_y[:,(rng<=max_rng)].ravel()]).transpose()
vals = cref_next[:,(rng<=max_rng)].ravel().data

cref_interp_flat = interpolate_to_points(pts, vals, xi, interp_type='linear', search_radius=hres)
cref_interp_next = cref_interp_flat.reshape(xm.shape)
cref_smooth_next = filters.gaussian(cref_interp_next, sigma=0.8)



#%% Find QLCS objects in previous and next volumes

qlcs_labels_prev, qlcs_regions_prev = find_qlcs_objects(cref_smooth_prev, hres, merge_distance=6)
qlcs_labels_next, qlcs_regions_next = find_qlcs_objects(cref_smooth_next, hres, merge_distance=6)



#%% Get centroids from analysis scan

# ib1,jb1,ib2,jb2 = qlcs_region.bbox
# bbox = [[xm[ib1,jb1], ym[ib1,jb1]], [xm[ib1,jb2-1], ym[ib1,jb2-1]], [xm[ib2-1,jb1], ym[ib2-1,jb1]], [xm[ib2-1,jb2-1], ym[ib2-1,jb2-1]]]

i_centroid,j_centroid = qlcs_region.centroid
ic1,jc1,ic2,jc2 = np.floor(i_centroid), np.floor(j_centroid), np.ceil(i_centroid), np.ceil(j_centroid)
i1,j1,i2,j2 = int(ic1), int(jc1), int(ic2), int(jc2)
x_centroid = (1 - abs(j_centroid-jc1))*xm[i1,j1] + (1 - abs(j_centroid-jc2))*xm[i1,j2]
y_centroid = (1 - abs(i_centroid-ic1))*ym[i1,j1] + (1 - abs(i_centroid-ic2))*ym[i2,j1]
centroid = [x_centroid, y_centroid]


#%% Get centroids from previous and next scan

qlcs_labels_prev_old = qlcs_labels_prev
qlcs_labels_next_old = qlcs_labels_next


qlcs_obj_prev = np.zeros(shape=qlcs_labels_prev.shape, dtype=int)
qlcs_obj_next = np.zeros(shape=qlcs_labels_next.shape, dtype=int)


if len(qlcs_regions_prev) > 1:
    xc = np.zeros((len(qlcs_regions_prev),))
    yc = np.zeros((len(qlcs_regions_prev),))
    centroid_dist2 = np.zeros((len(qlcs_regions_prev),))
    regs = np.unique(qlcs_labels_prev[(qlcs_labels_prev>0)])
    for n in range(len(centroid_dist2)):
        i_centroid,j_centroid = qlcs_regions_prev[n].centroid
        ic1,jc1,ic2,jc2 = np.floor(i_centroid), np.floor(j_centroid), np.ceil(i_centroid), np.ceil(j_centroid)
        i1,j1,i2,j2 = int(ic1), int(jc1), int(ic2), int(jc2)
        xc[n] = (1 - abs(j_centroid-jc1))*xm[i1,j1] + (1 - abs(j_centroid-jc2))*xm[i1,j2]
        yc[n] = (1 - abs(i_centroid-ic1))*ym[i1,j1] + (1 - abs(i_centroid-ic2))*ym[i2,j1]
        
        centroid_dist2[n] = (xc[n] - x_centroid)**2 + (yc[n] - y_centroid)**2
        
    n_closest = np.argmin(centroid_dist2)
    
    qlcs_obj_prev[(qlcs_labels_prev == n_closest+1)] = 1
    # qlcs_obj_prev[(qlcs_labels_prev == regs[n_closest])] = 1
    qlcs_region_prev = qlcs_regions_prev[n_closest]
    x_centroid_prev = xc[n_closest]
    y_centroid_prev = yc[n_closest]
else:
    i_centroid,j_centroid = qlcs_regions_prev[0].centroid
    ic1,jc1,ic2,jc2 = np.floor(i_centroid), np.floor(j_centroid), np.ceil(i_centroid), np.ceil(j_centroid)
    i1,j1,i2,j2 = int(ic1), int(jc1), int(ic2), int(jc2)
    x_centroid_prev = (1 - abs(j_centroid-jc1))*xm[i1,j1] + (1 - abs(j_centroid-jc2))*xm[i1,j2]
    y_centroid_prev = (1 - abs(i_centroid-ic1))*ym[i1,j1] + (1 - abs(i_centroid-ic2))*ym[i2,j1]
    
    qlcs_obj_prev = qlcs_labels_prev
    


if len(qlcs_regions_next) > 1:
    xc = np.zeros((len(qlcs_regions_next),))
    yc = np.zeros((len(qlcs_regions_next),))
    centroid_dist2 = np.zeros((len(qlcs_regions_next),))
    regs = np.unique(qlcs_labels_next[(qlcs_labels_next>0)])
    for n in range(len(centroid_dist2)):
        i_centroid,j_centroid = qlcs_regions_next[n].centroid
        ic1,jc1,ic2,jc2 = np.floor(i_centroid), np.floor(j_centroid), np.ceil(i_centroid), np.ceil(j_centroid)
        i1,j1,i2,j2 = int(ic1), int(jc1), int(ic2), int(jc2)
        xc[n] = (1 - abs(j_centroid-jc1))*xm[i1,j1] + (1 - abs(j_centroid-jc2))*xm[i1,j2]
        yc[n] = (1 - abs(i_centroid-ic1))*ym[i1,j1] + (1 - abs(i_centroid-ic2))*ym[i2,j1]
        
        centroid_dist2[n] = (xc[n] - x_centroid)**2 + (yc[n] - y_centroid)**2
        
    n_closest = np.argmin(centroid_dist2)
    
    qlcs_obj_next[(qlcs_labels_next == n_closest+1)] = 1
    # qlcs_obj_next[(qlcs_labels_next == regs[n_closest])] = 1
    qlcs_region_next = qlcs_regions_next[n_closest]
    x_centroid_next = xc[n_closest]
    y_centroid_next = yc[n_closest]
else:
    i_centroid,j_centroid = qlcs_regions_next[0].centroid
    ic1,jc1,ic2,jc2 = np.floor(i_centroid), np.floor(j_centroid), np.ceil(i_centroid), np.ceil(j_centroid)
    i1,j1,i2,j2 = int(ic1), int(jc1), int(ic2), int(jc2)
    x_centroid_next = (1 - abs(j_centroid-jc1))*xm[i1,j1] + (1 - abs(j_centroid-jc2))*xm[i1,j2]
    y_centroid_next = (1 - abs(i_centroid-ic1))*ym[i1,j1] + (1 - abs(i_centroid-ic2))*ym[i2,j1]
    
    qlcs_obj_next = qlcs_labels_next


#%% Estimate storm motion from centroid displacement

from datetime import datetime

time = cradar.time['mean']
time_prev = cradar_prev.time['mean']
time_next = cradar_next.time['mean']


dt_prev = (time - time_prev).total_seconds()
dt_next = (time_next - time).total_seconds()


u_prev = (x_centroid - x_centroid_prev)*1000 / dt_prev
v_prev = (y_centroid - y_centroid_prev)*1000 / dt_prev

u_next = (x_centroid_next - x_centroid)*1000 / dt_next
v_next = (y_centroid_next - y_centroid)*1000 / dt_next

u_sm = np.mean([u_prev, u_next])
v_sm = np.mean([v_prev, v_next])




#%% Leading line identification



























#%%


fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, qlcs_labels_prev, 'dbz', ax, datalims=[0,np.max(qlcs_labels_prev)], cmap='HomeyerRainbow')
ax.scatter(x_centroid_prev, y_centroid_prev, marker='.', s=30, c='k')
ax.set_title(f"Final QLCS objects\n {filenames[ind-2][101:]} ")
plt.show()

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, qlcs_labels, 'dbz', ax, datalims=[0,np.max(qlcs_labels)], cmap='HomeyerRainbow')
ax.scatter(x_centroid, y_centroid, marker='.', s=30, c='k')
ax.set_title(f"Final QLCS objects\n {filename[101:]} ")
plt.show()

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, qlcs_labels_next, 'dbz', ax, datalims=[0,np.max(qlcs_labels_next)], cmap='HomeyerRainbow')
ax.scatter(x_centroid_next, y_centroid_next, marker='.', s=30, c='k')
ax.set_title(f"Final QLCS objects\n {filenames[ind+2][101:]} ")
plt.show()




fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, qlcs_obj_prev, 'dbz', ax, datalims=[0,1], cmap='HomeyerRainbow')
ax.scatter(x_centroid_prev, y_centroid_prev, marker='.', s=30, c='k')
ax.set_title(f"Final QLCS objects\n {filenames[ind-2][101:]} ")
plt.show()

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, qlcs_obj, 'dbz', ax, datalims=[0,1], cmap='HomeyerRainbow')
ax.scatter(x_centroid, y_centroid, marker='.', s=30, c='k')
ax.set_title(f"Final QLCS objects\n {filename[101:]} ")
plt.show()

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, qlcs_obj_next, 'dbz', ax, datalims=[0,1], cmap='HomeyerRainbow')
ax.scatter(x_centroid_next, y_centroid_next, marker='.', s=30, c='k')
ax.set_title(f"Final QLCS objects\n {filenames[ind+2][101:]} ")
plt.show()




fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_smooth_prev, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
ax.set_title('Gridded CRef with Gaussian smoothing, prev')
plt.show()

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
ax.set_title('Gridded CRef with Gaussian smoothing')
plt.show()

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_smooth_next, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
ax.set_title('Gridded CRef with Gaussian smoothing, next')
plt.show()



























#%% QLCS object identification development!

# cref_bin = np.zeros(shape=cref.shape, dtype=bool)
# cref_bin[(cref>41)] = True
# cref_bin_clean = morphology.remove_small_objects(cref_bin, min_size=864, connectivity=1)
cref_bin = np.zeros(shape=cref_smooth.shape, dtype=bool)
cref_bin[(cref_smooth>40)] = True
cref_bin_clean = morphology.remove_small_objects(cref_bin, min_size=6, connectivity=1)


cref_labeled = measure.label(cref_bin_clean, connectivity=2)

regions = measure.regionprops(cref_labeled)




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



fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_bin, 'dbz', ax, datalims=[0,1], cmap='Grays')
ax.set_title('Binary CRef > 40 dBZ')
plt.show()


fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_bin_clean, 'dbz', ax, datalims=[0,1], cmap='Grays')
ax.set_title('Cleaned binary CRef > 40 dBZ')
plt.show()


fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_labeled, 'dbz', ax, datalims=[0,len(regions)], cmap='HomeyerRainbow')
# ax.scatter([xm[ix1,iy1], xm[ix1,iy2], xm[ix2,iy1], xm[ix2,iy2]], [ym[ix1,iy1], ym[ix1,iy2], ym[ix2,iy1], ym[ix2,iy2]], marker='.', s=30, c='k')
ax.set_title('CRef objects')
plt.show()




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
plot_cfill(xm, ym, cref_labeled_filtered, 'dbz', ax, datalims=[0,len(regions)], cmap='HomeyerRainbow')
ax.set_title('Filtered CRef objects with max CRef > 45 dBZ')
plt.show()



cref_bin_filtered = np.zeros(shape=cref_smooth.shape, dtype=bool)
cref_bin_filtered[(cref_labeled_filtered>0)] = True

footprint = morphology.footprint_rectangle((9,9))
cref_bin_dilated = morphology.binary_dilation(cref_bin_filtered, footprint=footprint)

cref_merged = measure.label(cref_bin_dilated, connectivity=2)

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_bin_dilated, 'dbz', ax, datalims=[0,1], cmap='Grays')
ax.set_title('Dilated binary filtered CRef objects')
plt.show()


fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_merged, 'dbz', ax, datalims=[0,len(regions)], cmap='HomeyerRainbow')
ax.set_title('Merged filtered CRef objects')
plt.show()



# new_regions_filtered = []
# for i in regions_filtered:
#     new_regions_filtered.append( np.unique(cref_merged[(cref_labeled_filtered==i)])[0] )

# cref_merged_undilated = np.zeros(shape=cref_merged.shape, dtype=int)
# for i in range(len(regions_filtered)):
#     cref_merged_undilated[(cref_labeled_filtered==regions_filtered[i])] = new_regions_filtered[i]

# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, cref_merged_undilated, 'dbz', ax, datalims=[0,len(regions_filtered)], cmap='HomeyerRainbow')
# ax.set_title('Undilated merged filtered CRef objects')
# plt.show()

cref_merged[(cref_bin_filtered==0)] = 0

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_merged, 'dbz', ax, datalims=[0,len(regions)], cmap='HomeyerRainbow')
ax.set_title('Undilated merged filtered CRef objects')
plt.show()



regions_merged = measure.regionprops(cref_merged)



regions_fin = []
cref_merged_filtered = np.zeros(shape=cref_merged.shape, dtype=int)

for i in range(len(regions_merged)):
    area = regions_merged[i].area
    axis_major = regions_merged[i].axis_major_length
    ecc = regions_merged[i].eccentricity
    
    print(f"Region {i+1} major axis= {axis_major:.1f} pixels ({axis_major*3:.1f} km) , ecc={ecc:.2f}")
    
    
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
        # cref_merged_filtered[(cref_merged==i+1)] = i+1
        cref_merged_filtered[(cref_merged==i+1)] = i+1
        # cref_merged_filtered[(cref_bin_filtered==0)] = 0
        regions_fin.append(i+1)
        


### YESSSS IT WORKS

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_merged_filtered, 'dbz', ax, datalims=[0,np.max(regions_fin)], cmap='HomeyerRainbow')
ax.set_title('Final QLCS objects\n filtered for max dBZ, merged, and filtered for length and eccentricity ')
plt.show()







