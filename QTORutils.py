# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 11:41:04 2026

@author: mschne28

Modules and functions for calculating QTor.
"""

####################
### Load modules ###
####################

import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import pyart #need an earlier version of xarray -> 0.20.2 or earlier
import pickle
# import xarray as xr
from skimage import measure,morphology,filters
# import sklearn
# from glob import glob
import os
from os.path import exists
from matplotlib.ticker import MultipleLocator

#%%

# Define colormaps for common cm1 variables
cmaps = {
    'u':       {'cm': 'balance',        'label': "u (m s$^{-1}$)"},
    'v':       {'cm': 'balance',        'label': "v (m s$^{-1}$)"},
    'w':       {'cm': 'balance',        'label': "w (m s$^{-1}$)"},
    'wspd':    {'cm': 'HomeyerRainbow', 'label': "Wind speed (m s$^{-1}$)"},
    'th':      {'cm': 'HomeyerRainbow', 'label': "\u03B8 (K)"},
    'thpert':  {'cm': 'balance',        'label': "\u03B8' (K)"},
    'thr':     {'cm': 'HomeyerRainbow', 'label': "\u03B8\u1D68 (K)"},
    'thrpert': {'cm': 'balance',        'label': "\u03B8'\u1D68 (K)"},
    'qv':      {'cm': 'YlGnBu',         'label': "q$_v$ (g kg$^{-1}$)"},
    'qvpert':  {'cm': 'balance',        'label': "q'$_v$ (g kg$^{-1}$)"},
    'qc':      {'cm': 'BuPu',           'label': "q$_c$ (g kg$^{-1}$)"},
    'qr':      {'cm': 'BuPu',           'label': "q$_r$ (g kg$^{-1}$)"},
    'qi':      {'cm': 'BuPu',           'label': "q$_i$ (g kg$^{-1}$)"},
    'prs':     {'cm': 'HomeyerRainbow', 'label': "p (hPa)"},
    'prspert': {'cm': 'balance',        'label': "p' (hPa)"},
    'pi':      {'cm': 'HomeyerRainbow', 'label': "\u03C0 (nondimensional)"},
    'pipert':  {'cm': 'balance',        'label': "\u03C0' (nondimensional)"},
    'rho':     {'cm': 'HomeyerRainbow', 'label': "\u03C1 (kg m$^{-3}$)"},
    'xvort':   {'cm': 'balance',        'label': "\u03BE (s$^{-1}$)"},
    'yvort':   {'cm': 'balance',        'label': "\u03B7 (s$^{-1}$)"},
    'zvort':   {'cm': 'balance',        'label': "\u03B6 (s$^{-1}$)"},
    'hvort':   {'cm': 'HomeyerRainbow', 'label': "\u03c9$_H$ (s$^{-1}$)"},
    'vort':    {'cm': 'HomeyerRainbow', 'label': "\u03c9 (s$^{-1}$)"},
    'OW':      {'cm': 'balance',        'label': "OW (s$^{-2}$)"},
    'divh':    {'cm': 'balance',        'label': "\u25BD$_H$u (s$^{-1}$)"},
    'dbz':     {'cm': 'HomeyerRainbow', 'label': "$Z_H$ (dBZ)"},
    'pgf':     {'cm': 'balance',        'label': "PPGA (m s$^{-2}$)"},
    'cape':    {'cm': 'HomeyerRainbow', 'label': "CAPE (J kg$^{-1}$)"},
    'srh':     {'cm': 'HomeyerRainbow', 'label': "SRH (m$^{2}$ s$^{-2}$)"},
    'z':       {'cm': 'HomeyerRainbow', 'label': "Height (m)"},
    'scp':     {'cm': 'HomeyerRainbow', 'label': "SCP"},
    'stp':     {'cm': 'HomeyerRainbow', 'label': "STP"},
    'uh':      {'cm': 'HomeyerRainbow', 'label': "UH (m$^{2}$ s$^{-2}$)"},
    'del2thp': {'cm': 'balance',        'label': "\u25BD$^2$\u03B8' (K m$^{-2}$)"}
}



#%% Data processing and analysis functions


def cressman_interpolation(obs_lats, obs_lons, data, grid_lats, grid_lons, radius):
    '''
    obs_lats  : 2D array-like of observation lats (y) - shape (M,N)
    obs_lons  : 2D array-like of observation lons (x) - shape (M,N)
    data      : 2D array-like of observation values - shape (M,N)
    grid_lats : 2D array-like of grid lats - shape (P,Q)
    grid_lons : 2D array-like of grid lons - shape (P,Q)
    radius    : Radius of influence
    '''
    
    # Convert inputs to numpy arrays
    obs_lats = np.asarray(obs_lats, dtype=np.float32)
    obs_lons = np.asarray(obs_lons, dtype=np.float32)
    data = np.asarray(data, dtype=np.float32)
    grid_lats = np.asarray(grid_lats, dtype=np.float32)
    grid_lons = np.asarray(grid_lons, dtype=np.float32)
    
    # Flatten grid for vectorized computation
    obs_lat_flat = obs_lats.ravel()
    obs_lon_flat = obs_lons.ravel()
    data_flat = data.ravel()
    grid_lat_flat = grid_lats.ravel()
    grid_lon_flat = grid_lons.ravel()
    
    # Compute squared distances between each grid point and each observation
    # Broadcasting: (G,1) - (1,O) → (G,O)
    dlat2 = (grid_lat_flat[:,None] - obs_lat_flat[None,:])**2
    dlon2 = (grid_lon_flat[:,None] - obs_lon_flat[None,:])**2
    dist2 = dlat2 + dlon2
    
    # Apply Cressman weights only where dist2 < radius^2
    mask = dist2 < radius**2
    weights = np.zeros_like(dist2)
    weights[mask] = (radius**2 - dist2[mask]) / (radius**2 + dist2[mask])
    
    # Weighted sum and normalization
    weighted_sum = np.sum(weights*data_flat[None,:], axis=1)
    weight_total = np.sum(weights, axis=1)
    
    # Avoid division by 0
    with np.errstate(divide='ignore', invalid='ignore'):
        analysis_flat = np.where(weight_total>0, weighted_sum/weight_total, np.nan)
    
    # Reshape data back to grid shape
    analysis = analysis_flat.reshape(grid_lats.shape)
    
    return analysis
    
    
    


# Retrieve QLCS objects

def find_qlcs_objects(cref, hres, min_cref=40, max_cref=45, merge_distance=12, min_length1=100, min_length2=150, min_ecc1=0.85, min_ecc2=0.74, min_area=54):
    '''
    QLCS object identification following Britt et al. 2024 and 2026.
    
    cref : Composite reflectivity interpolated onto 2D Cartesian grid.
    hres : Horizontal grid resolution.
    min_cref : First reflectivity threshold for initial storm object ID. Default is 40 dBZ.
    max_cref : Second reflectivity threshold for filtering storm objects. Default is 45 dBZ.
    merge_distance : Distance threshold for merging storm objects. Default is 12 km.
    min_length1 : First length threshold for QLCS storm objects. Default is 100 km.
    min_length2 : Second length threshold for QLCS storm objects. Default is 150 km.
    min_ecc1 : First eccentricity threshold for QLCS storm objects, corresponding to objects with length >length1 and <length2. Default is 0.85.
    min_ecc2 : Second eccentricity threshold for QLCS storm objects, corresponding to objects with length >length2. Default is 0.74.
    min_area : Area threshold for storm objects. Default is 54 km^2.
    '''
    
    # Binarize smoothed cref according to min_cref
    cref_bin = np.zeros(shape=cref.shape, dtype=bool)
    cref_bin[(cref > min_cref)] = True
    
    # Clean up small areas less than min_area
    min_size = min_area / (hres**2) #convert to number of pixels
    cref_bin_clean = morphology.remove_small_objects(cref_bin, min_size=min_size, connectivity=1)
    
    # Generate object labels for each binary object
    obj_labels = measure.label(cref_bin_clean, connectivity=2)
    
    # Generate region properties for each object
    regions = measure.regionprops(obj_labels)
    
    
    # Filter initial storm objects with max reflectivity below second threshold
    regions_filtered = []
    obj_labels_filtered = np.zeros(shape=obj_labels.shape, dtype=int)

    for i in range(len(regions)):
        area = regions[i].area
        axis_major = regions[i].axis_major_length
        coords = regions[i].coords
        bbox = regions[i].bbox
        ecc = regions[i].eccentricity
        
        if np.nanmax(cref[(obj_labels==i+1)]) > max_cref:
            maxz_met = True
            obj_labels_filtered[(obj_labels==i+1)] = i+1
            regions_filtered.append(i+1)
        else:
            maxz_met = False
    
    
    # Rebinarize filtered object labels
    cref_bin_filtered = np.zeros(shape=cref.shape, dtype=bool)
    cref_bin_filtered[(obj_labels_filtered>0)] = True
    
    # Dilate binarized objects to merge objects within merge_distance of each other
    merge_len_oneway = merge_distance / hres #convert to number of pixels
    if np.mod(merge_distance, hres) != 0:
        merge_len_oneway = np.round(merge_len_oneway)
    merge_len = int(2*merge_len_oneway + 1)
    footprint = morphology.footprint_rectangle((merge_len, merge_len))
    cref_bin_dilated = morphology.binary_dilation(cref_bin_filtered, footprint=footprint)
    
    # Retrieve new labels for dilated filtered binary field, then undilate back to original object shapes
    obj_labels_merged = measure.label(cref_bin_dilated, connectivity=2)
    obj_labels_merged[(cref_bin_filtered==0)] = 0
    
    # Generate new objects from the merged filtered binary field
    regions_merged = measure.regionprops(obj_labels_merged)
    
    
    # Filter merged regions 
    regions_final = []
    qlcs_labels_final = np.zeros(shape=obj_labels_merged.shape, dtype=int)
    qlcs_regions_final = []
    
    n = 0
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
            n = n+1
            # isQLCS = True
            qlcs_labels_final[(obj_labels_merged==i+1)] = n
            qlcs_labels_final[(cref_bin_filtered==0)] = 0
            regions_final.append(i+1)
            
            qlcs_regions_final.append(regions_merged[i])
    
    return qlcs_labels_final, qlcs_regions_final
    
    
    
    
    


#%% Miscellaneous other functions


# Convert lat/lon coordinates to x/y distances relative to an origin point (in km)
def latlon2xy(lat, lon, lat_o, lon_o):
    # lat, lon:     1-D vectors of lat/lon in decimal degrees N/deg E
    # lat_o, lon_o: lat/lon of origin in decimal degrees N/deg E
    
    r_earth = 6378.1 # km
    
    thy = lat_o*np.pi/180 # convert to radians
    thz = -lon_o*np.pi/180
    
    # transform matrices
    Ry = [[np.cos(thy),  0,  np.sin(thy)],
          [0,            1,  0],
          [-np.sin(thy), 0,  np.cos(thy)]]
    
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



# Save data to new or existing pickle file
def save_to_pickle(data, pkl_fname, new_pkl=False):
    # data: dict of variables to save
    # pkl_fname: filename to save data to (includes path)
    # new_pkl: True or False (if True, will overwrite any existing file)
    if (not exists(pkl_fname)) | new_pkl:
        dbfile = open(pkl_fname, 'wb')
        pickle.dump(data, dbfile)
        dbfile.close()
    elif exists(pkl_fname):
        dbfile = open(pkl_fname, 'rb')
        save_data = pickle.load(dbfile)
        dbfile.close()
        
        save_data.update(data)
        dbfile = open(pkl_fname, 'wb')
        pickle.dump(save_data, dbfile)
        dbfile.close()



#%% Plotting functions

# Wrapper function for contourf
def plot_contourf(x, y, data, field, ax, levels=None, datalims=None, xlims=None, ylims=None,
                  cmap=None, cbar=True, cbfs=None, cbticks=None, **kwargs):
    if cmap is None:
        cm, cb_label = cmaps[field]['cm'], cmaps[field]['label']
    else:
        cm, cb_label = cmap, cmaps[field]['label']
    
    if levels is None:
        levs = None
    else:
        levs = levels
    
    if datalims is None:
        datamin = None
        datamax = None
    else:
        datamin = datalims[0]
        datamax = datalims[1]
    
    c = ax.contourf(x, y, data, levels=levs, vmin=datamin, vmax=datamax, cmap=cm, antialiased=True, **kwargs)
    # ax.contour(x, y, data, levels=levs, vmin=datamin, vmax=datamax, cmap=cm, antialiased=True, **kwargs)
    # ax.contourf(x, y, data, levels=levs, vmin=datamin, vmax=datamax, cmap=cm, antialiased=True, **kwargs)
    # ax.contourf(x, y, data, levels=levs, vmin=datamin, vmax=datamax, cmap=cm, antialiased=True, **kwargs)
    c.set_edgecolor('face')
    
    if cbar:
        cb = plt.colorbar(c, ax=ax, extend='both')
        # cb.set_label(cb_label)
        if np.max(np.abs(datalims)) < 0.1:
            cb.formatter.set_powerlimits((0,0))
        if cbfs is None:
            cb.set_label(cb_label)
        else:
            cb.set_label(cb_label, fontsize=cbfs)
        if cbticks is not None:
            cb.set_ticks(cbticks)
    
    if xlims is not None:
        ax.set_xlim(xlims[0], xlims[1])
    if ylims is not None:
        ax.set_ylim(ylims[0], ylims[1])
    
    return c



# Wrapper function for pcolormesh
def plot_cfill(x, y, data, field, ax, datalims=None, xlims=None, ylims=None,
               cmap=None, cbar=True, cbfs=None, cbticks=None, **kwargs):
    if cmap is None:
        cm, cb_label = cmaps[field]['cm'], cmaps[field]['label']
    else:
        cm, cb_label = cmap, cmaps[field]['label']
    
    if datalims is None:
        datamin = None
        datamax = None
    else:
        datamin = datalims[0]
        datamax = datalims[1]
    
    # Create the plot
    c = ax.pcolormesh(x, y, data, vmin=datamin, vmax=datamax, cmap=cm, **kwargs)

    # Format the colorbar
    # c.cmap.set_bad('grey', 1.0)
    if cbar:
        cb = plt.colorbar(c, ax=ax, extend='both')
        cb.set_label(cb_label)
        if np.max(np.abs(datalims)) < 0.1:
            cb.formatter.set_powerlimits((0,0))
        if cbfs is None:
            cb.set_label(cb_label)
        else:
            cb.set_label(cb_label, fontsize=cbfs)
        if cbticks is not None:
            cb.set_ticks(cbticks)
    
    if xlims is not None:
        ax.set_xlim(xlims[0], xlims[1])
    if ylims is not None:
        ax.set_ylim(ylims[0], ylims[1])
    
    return c

















