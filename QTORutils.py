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
import xarray as xr
# import sklearn
from glob import glob
import os
from os.path import exists
from matplotlib.ticker import MultipleLocator
from scipy.ndimage import gaussian_filter
from skimage import measure,morphology,filters
from metpy.interpolate import interpolate_to_points
from datetime import datetime

from scipy.ndimage import gaussian_filter1d
import math

import shapely
from shapely.geometry import Polygon,Point,LineString,MultiLineString
from shapely.plotting import plot_polygon
# from centerline.geometry import Centerline
from shapely.ops import linemerge,substring
import networkx as nx

from label_centerlines import get_centerline

from scipy.spatial import KDTree


# import geojson
# import geopandas as gpd



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


# raster?
def calc_QTor():
    return



# Single point
def get_local_tortuosity():
    return



# Single point
def get_line_parallel_shear(shear, theta_local, theta_norm):
    return



# Single point
def get_line_normal_shear(shear, theta_local, theta_norm):
    return



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



# Get ERA5 shear components (raster)
def get_shear(z, orog, u, v, u10, v10, depth, bottom=None):
    if bottom is None:
        bottom = 10
        top = depth
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
    
    u_shear = u_top - u10
    v_shear = v_top - v10
    
    # mask = ((z_interp>bottom) | np.isclose(z_interp,bottom)) & ((z_interp<top) | np.isclose(z_interp,top))
    # # z_masked = np.ma.masked_array(z_interp, ~mask)
    # u_layer = np.ma.masked_array(u_interp, ~mask)
    # v_layer = np.ma.masked_array(v_interp, ~mask)
    # # u_layer = u_interp[(mask)]
    # # v_layer = v_interp[(mask)]
    
    # u_shear = u_layer[-1,:,:] - u_layer[0,:,:]
    # v_shear = v_layer[-1,:,:] - v_layer[0,:,:]
    
    return (u_shear, v_shear)



# Get inflow polygon
def get_inflow_polygon(ll_points, storm_motion):
    # ll_points: array-like, column 1: x coords, column 2: y coords
    # storm_motion: array-like, [u_velocity, v_velocity] of storm object
    
    x_ll = ll_points[:,0]
    y_ll = ll_points[:,1]
    u_storm = storm_motion[0]
    v_storm = storm_motion[1]
    
    # Project leading line forward
    theta_local = np.zeros((len(ll_points),), dtype=float)
    theta_norm = np.zeros((len(ll_points),), dtype=float)
    x_proj = np.zeros((len(ll_points),), dtype=float)
    y_proj = np.zeros((len(ll_points),), dtype=float)
    dist = 40 #km, near-field distance

    for i in range(len(ll_points)):
        if i == 0:
            theta_local[i] = np.arctan2(y_ll[i+1]-y_ll[i], x_ll[i+1]-x_ll[i])
        elif i == len(ll_points)-1:
            theta_local[i] = np.arctan2(y_ll[i]-y_ll[i-1], x_ll[i]-x_ll[i-1])
        else:
            dist1 = distance([x_ll[i],y_ll[i]], [x_ll[i-1],y_ll[i-1]])
            dist2 = distance([x_ll[i],y_ll[i]], [x_ll[i+1],y_ll[i+1]])
            theta1 = np.arctan2(y_ll[i]-y_ll[i-1], x_ll[i]-x_ll[i-1])
            theta2 = np.arctan2(y_ll[i+1]-y_ll[i], x_ll[i+1]-x_ll[i])
            theta_local[i] = dist2/(dist1+dist2)*theta1 + dist1/(dist1+dist2)*theta2

    # Determine whether to add or subtract 90 deg from local angle to get normal angle
    for i in range(len(theta_local)):
        if abs(u_storm) > abs(v_storm):
            if u_storm > 0:
                theta_norm[i] = theta_local[i] + np.pi/2
            elif u_storm < 0:
                theta_norm[i] = theta_local[i] - np.pi/2
        elif abs(v_storm) > abs(u_storm):
            if v_storm > 0:
                theta_norm[i] = theta_local[i] + np.pi/2
            elif v_storm < 0:
                theta_norm[i] = theta_local[i] - np.pi/2
        
        x_proj[i] = x_ll[i] + dist*np.cos(theta_norm[i])
        y_proj[i] = y_ll[i] + dist*np.sin(theta_norm[i])
    
    # QC polygon - make sure all projection points are ahead of leading line
    qc_points = ll_points.tolist()
    qc_points.append([-300,-300])
    qc_polygon = Polygon(qc_points)

    theta_norm_qc = np.zeros((len(ll_points),))
    x_proj_qc = np.zeros((len(ll_points),))
    y_proj_qc = np.zeros((len(ll_points),))

    for i in range(len(ll_points)):
        point = Point(x_proj[i], y_proj[i])
        is_inside = qc_polygon.contains(point)
        if is_inside:
            if theta_norm[i] > 0:
                theta_norm_qc[i] = theta_norm[i] - np.pi
            else:
                theta_norm_qc[i] = theta_norm[i] + np.pi
            x_proj_qc = x_ll[i] + dist*np.cos(theta_norm_qc[i])
            y_proj_qc = y_ll[i] + dist*np.sin(theta_norm_qc[i])
        else:
            theta_norm_qc[i] = theta_norm[i]
            x_proj_qc[i] = x_proj[i]
            y_proj_qc[i] = y_proj[i]

    leading_points = ll_points.tolist()
    proj_points = [[x_proj_qc[i], y_proj_qc[i]] for i in range(len(ll_points))][::-1]

    inflow_polygon = shapely.buffer(Polygon(leading_points+proj_points), 0)
    
    return inflow_polygon,theta_norm,theta_local



def distance(point1, point2):
    return math.hypot(point2[0] - point1[0], point2[1] - point1[1])



# Find leading line
def get_leading_line(qlcs_obj, storm_motion, grid_x, grid_y, grid_lat, grid_lon):
    # qlcs_obj: array-like, bool of object
    # storm_motion: array-like, [u_velocity, v_velocity] of storm object
    # grid_x, grid_y: x/y meshgrids
    # grid_lat, grid_lon: lat/lon meshgrids
    
    # Initial leading line first guess
    u_storm = storm_motion[0]
    v_storm = storm_motion[1]
    
    # E-W search
    x_ll_zonal = []
    y_ll_zonal = []
    lon_ll_zonal = []
    lat_ll_zonal = []
    for jdy in range(grid_x.shape[0]):
        if np.any(qlcs_obj[jdy,:] > 0):
            if u_storm > 0:
                idx = np.where(qlcs_obj[jdy,:]>0)[0][-1]
            elif u_storm < 0:
                idx = np.where(qlcs_obj[jdy,:]>0)[0][0]
            x_ll_zonal.append(grid_x[jdy,idx])
            y_ll_zonal.append(grid_y[jdy,idx])
            lon_ll_zonal.append(grid_lon[jdy,idx])
            lat_ll_zonal.append(grid_lat[jdy,idx])
            

    # N-S search
    x_ll_merid = []
    y_ll_merid = []
    lon_ll_merid = []
    lat_ll_merid = []
    for idx in range(grid_y.shape[1]):
        if np.any(qlcs_obj[:,idx] > 0):
            if v_storm > 0:
                jdy = np.where(qlcs_obj[:,idx]>0)[0][-1]
            elif v_storm < 0:
                jdy = np.where(qlcs_obj[:,idx]>0)[0][0]
            x_ll_merid.append(grid_x[jdy,idx])
            y_ll_merid.append(grid_y[jdy,idx])
            lon_ll_merid.append(grid_lon[jdy,idx])
            lat_ll_merid.append(grid_lat[jdy,idx])
    
    ### QC ###
    
    # Find centerline
    qlcs_obj_smooth = filters.gaussian(qlcs_obj, sigma=5/3)
    thres = np.percentile(qlcs_obj_smooth[(qlcs_obj_smooth>0)], 60) #could just use 5e-20? idk if the limits are universal
    obj_smooth_bin = np.where(qlcs_obj_smooth>thres, 1, 0)

    contours = measure.find_contours(obj_smooth_bin.astype(bool), level=0.99, fully_connected='low')
    ind = np.asarray([len(contours[n]) for n in range(len(contours))])
    contour = contours[np.argmax(ind)]
    obj_contour_int = np.round(contour).astype(int)
    obj_contour = np.array([ [grid_x[j,i], grid_y[j,i]] for j,i in zip(obj_contour_int[:,0], obj_contour_int[:,1]) ])
    
    polygon = Polygon(obj_contour)
    centerline = get_centerline(polygon)
    trimmed_centerline = substring(centerline, 0.05*centerline.length, 0.95*centerline.length)
    x_cl,y_cl = trimmed_centerline.xy
    
    # QC leading line
    cl_points = np.column_stack( (x_cl, y_cl))
    zonal_points = np.column_stack( (x_ll_zonal, y_ll_zonal))
    merid_points = np.column_stack( (x_ll_merid, y_ll_merid))
    zonal_lonlat = np.column_stack( (lon_ll_zonal, lat_ll_zonal))
    merid_lonlat = np.column_stack( (lon_ll_merid, lat_ll_merid))
    
    tree_cl = KDTree(cl_points)
    dist_zonal,inds_zonal = tree_cl.query(zonal_points) # Find nearest centerline point in tree_cl to each zonal point
    dist_merid,inds_merid = tree_cl.query(merid_points) # Find nearest centerline point in tree_cl to each meridional point

    cl_points_matched_zonal = np.array([cl_points[inds_zonal[i]] for i in range(len(inds_zonal))])
    cl_points_matched_merid = np.array([cl_points[inds_merid[i]] for i in range(len(inds_merid))])

    vector_zonal = zonal_points - cl_points_matched_zonal
    vector_merid = merid_points - cl_points_matched_merid
    
    if abs(u_storm) > abs(v_storm):
        zonal_points_filtered = zonal_points[(vector_zonal[:,0]/u_storm > 0),:]
        merid_points_filtered = merid_points[(vector_merid[:,0]/u_storm > 0),:]
        zonal_lonlat_filtered = zonal_lonlat[(vector_zonal[:,0]/u_storm > 0),:]
        merid_lonlat_filtered = merid_lonlat[(vector_merid[:,0]/u_storm > 0),:]
    elif abs(v_storm) > abs(u_storm):
        zonal_points_filtered = zonal_points[(vector_zonal[:,1]/v_storm > 0),:]
        merid_points_filtered = merid_points[(vector_merid[:,1]/v_storm > 0),:]
        zonal_lonlat_filtered = zonal_lonlat[(vector_zonal[:,1]/v_storm > 0),:]
        merid_lonlat_filtered = merid_lonlat[(vector_merid[:,1]/v_storm > 0),:]
    
    
    # Merge zonal and meridional points
    ll_points_cat = []
    ll_lonlat_cat = []
    seen = set()
    for item,item2 in zip(zonal_points_filtered.tolist() + merid_points_filtered.tolist(), zonal_lonlat_filtered.tolist() + merid_lonlat_filtered.tolist()):
        # if tuple(ll_points_cat[i]) not in seen:
        if tuple(item) not in seen:
            seen.add(tuple(item))
            ll_points_cat.append(item)
            ll_lonlat_cat.append(item2)

    ll_points_merged = []
    ll_lonlat_merged = []
    # query_point = [-300, 300]
    points_tmp = [item for item in ll_points_cat]
    lonlat_tmp = [item for item in ll_lonlat_cat]
    # is_sorted = [False] * len(ll_points_cat)
    tree_ll = KDTree(np.asarray(ll_points_cat))
    dist,ind = tree_ll.query([-300,300])
    start_point = ll_points_cat[ind]
    points_tmp.remove(start_point)
    start_lonlat = ll_lonlat_cat[ind]
    lonlat_tmp.remove(start_lonlat)
    # query_point = start_point
    ll_points_merged.append(start_point)
    ll_lonlat_merged.append(start_lonlat)

    for i in range(len(ll_points_cat)-1):
        tree = KDTree(np.asarray(points_tmp))
        dist,ind = tree.query(start_point)
        # dist2,ind2 = tree.query(query_point)
        ll_points_merged.append(points_tmp[ind])
        points_tmp.remove(points_tmp[ind])
        ll_lonlat_merged.append(lonlat_tmp[ind])
        lonlat_tmp.remove(lonlat_tmp[ind])

    ll_points_sorted = [item for item in ll_points_merged]
    ll_lonlat_sorted = [item for item in ll_lonlat_merged]
    # fix messed up points
    for i in range(1, len(ll_points_cat)-1):
        d1 = distance(ll_points_sorted[i-1], ll_points_sorted[i])
        d2 = distance(ll_points_sorted[i], ll_points_sorted[i+1])
        d3 = distance(ll_points_sorted[i-1], ll_points_sorted[i+1])
        
        if (d1>d3) & (d2>d3):
            ll_points_sorted[i], ll_points_sorted[i+1] = ll_points_merged[i+1], ll_points_merged[i]
            ll_lonlat_sorted[i], ll_lonlat_sorted[i+1] = ll_lonlat_merged[i+1], ll_lonlat_merged[i]
    
    # Smooth leading line
    ll_points_sorted = np.asarray(ll_points_sorted, dtype=float)
    ll_points_smooth = gaussian_filter1d(ll_points_sorted, sigma=2/3, axis=0)
    
    ll_lonlat_sorted = np.asarray(ll_lonlat_sorted, dtype=float)
    ll_lonlat_smooth = gaussian_filter1d(ll_lonlat_sorted, sigma=2/3, axis=0)
    
    return ll_points_smooth, ll_lonlat_smooth




# Get storm object centroid motion
def get_storm_motion(centroid1, centroid2, time1, time2):
    # centroid1, centroid2:    array-like, [x_coordinate, y_coordinate] of centroids (in km)
    # time1, time2:    datetime objects - mean volume scan times
    
    dt = (time2 - time1).total_seconds()
    u_centroid = (centroid2[0] - centroid1[0])*1000 / dt
    v_centroid = (centroid2[1] - centroid1[1])*1000 / dt
    
    # if time2 > time1:
    #     dt = (time2 - time1).total_seconds()
    #     u_centroid = (centroid2[0] - centroid1[0])*1000 / dt
    #     v_centroid = (centroid2[1] - centroid1[1])*1000 / dt
    # else:
    #     dt = (time1 - time2).total_seconds()
    #     u_centroid = (centroid1[0] - centroid2[0])*1000 / dt
    #     v_centroid = (centroid1[1] - centroid2[1])*1000 / dt
    
    return u_centroid,v_centroid



# Get QLCS object centroid
def get_centroid(region, grid_x, grid_y, reference_centroid=None):
    # region:    RegionsProps object or list of RegionProps objects
    # grid_x, grid_y:    x/y meshgrids
    # reference_centroid: array-like, [x_coordinate, y_coordinate] of reference centroid
    if (len(region) > 1) and (reference_centroid is not None):
        xc = np.zeros((len(region),))
        yc = np.zeros((len(region),))
        x_centroid_ref = reference_centroid[0]
        y_centroid_ref = reference_centroid[1]
        centroid_dist2 = np.zeros((len(region),))
        for n in range(len(centroid_dist2)):
            j_centroid,i_centroid = region[n].centroid
            jc1,ic1,jc2,ic2 = np.floor(j_centroid), np.floor(i_centroid), np.ceil(j_centroid), np.ceil(i_centroid)
            j1,i1,j2,i2 = int(jc1), int(ic1), int(jc2), int(ic2)
            xc[n] = (1 - abs(i_centroid-ic1))*grid_x[j1,i1] + (1 - abs(i_centroid-ic2))*grid_x[j1,i2]
            yc[n] = (1 - abs(j_centroid-jc1))*grid_y[j1,i1] + (1 - abs(j_centroid-jc2))*grid_y[j2,i1]
            
            centroid_dist2[n] = (xc[n] - x_centroid_ref)**2 + (yc[n] - y_centroid_ref)**2
            
        n_closest = np.argmin(centroid_dist2)
        # qlcs_obj[(qlcs_labels == n_closest+1)] = 1
        x_centroid = xc[n_closest]
        y_centroid = yc[n_closest]
        
    elif (len(region) > 1) and (reference_centroid is None):
        print("More than one region --> Need a reference centroid")
        return
    
    else:
        if isinstance(region, list):
            j_centroid,i_centroid = region[0].centroid
        else:
            j_centroid,i_centroid = region.centroid
        jc1,ic1,jc2,ic2 = np.floor(j_centroid), np.floor(i_centroid), np.ceil(j_centroid), np.ceil(i_centroid)
        j1,i1,j2,i2 = int(jc1), int(ic1), int(jc2), int(ic2)
        x_centroid = (1 - abs(i_centroid-ic1))*grid_x[j1,i1] + (1 - abs(i_centroid-ic2))*grid_x[j1,i2]
        y_centroid = (1 - abs(j_centroid-jc1))*grid_y[j1,i1] + (1 - abs(j_centroid-jc2))*grid_y[j2,i1]
        # qlcs_obj = qlcs_labels
    
    centroid = [x_centroid, y_centroid]
    
    return centroid




# Retrieve QLCS objects
def get_qlcs_objects(cref, hres, min_cref=40, max_cref=45, merge_distance=12, min_length1=100, min_length2=150, min_ecc1=0.85, min_ecc2=0.74, min_area=54):
    '''
    QLCS object identification following Britt et al. 2024 and 2026.
    
    cref : Composite reflectivity interpolated onto 2D Cartesian grid (and smoothed if needed, ideally).
    hres : Horizontal grid resolution (assumed km).
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
    merge_len_oneway = 0.5*merge_distance / hres #convert to number of pixels
    if np.mod(0.5*merge_distance, hres) != 0:
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
        # area = regions_merged[i].area
        axis_major = regions_merged[i].axis_major_length
        ecc = regions_merged[i].eccentricity
        
        # print(f"Region {i+1} major axis= {axis_major:.1f} pixels ({axis_major*3:.1f} km) , ecc={ecc:.2f}")
        
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




# Cressman filtering interpolation - does not work with large datasets
def cressman_interpolation(obs_x, obs_y, data, grid_x, grid_y, radius):
    '''
    obs_x  : array-like of observation x points
    obs_y  : array-like of observation y points
    data   : array-like of observations
    grid_x : array-like of grid x points
    grid_y : array-like of grid y points
    radius : Radius of influence
    '''
    
    # Convert inputs to numpy arrays
    obs_x = np.asarray(obs_x, dtype=np.float32)
    obs_y = np.asarray(obs_y, dtype=np.float32)
    data = np.asarray(data, dtype=np.float32)
    grid_x = np.asarray(grid_x, dtype=np.float32)
    grid_y = np.asarray(grid_y, dtype=np.float32)
    
    # Flatten grid for vectorized computation
    obs_x_flat = obs_x.ravel()
    obs_y_flat = obs_y.ravel()
    data_flat = data.ravel()
    grid_x_flat = grid_x.ravel()
    grid_y_flat = grid_y.ravel()
    
    # Compute squared distances between each grid point and each observation
    # Broadcasting: (G,1) - (1,O) → (G,O)
    dx2 = (grid_x_flat[:,None] - obs_x_flat[None,:])**2
    dy2 = (grid_y_flat[:,None] - obs_y_flat[None,:])**2
    dist2 = dx2 + dy2
    
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
    analysis = analysis_flat.reshape(grid_x.shape)
    
    return analysis




# Cressman filtering interpolation - parallelized with Dask to work with large datasets 
def cressman_interpolation_dask(obs_x, obs_y, data, grid_x, grid_y, radius, chunk_size=5000):
    import dask.array as da
    from dask.diagnostics import ProgressBar
    
    obs_x = np.asarray(obs_x, dtype=np.float32)
    obs_y = np.asarray(obs_y, dtype=np.float32)
    data = np.asarray(data, dtype=np.float32)
    grid_x = np.asarray(grid_x, dtype=np.float32)
    grid_y = np.asarray(grid_y, dtype=np.float32)
    
    # Flatten obs and grid for vectorized computation
    obs_points = np.column_stack( (obs_x.ravel(), obs_y.ravel()) )
    obs_values = data.ravel()
    grid_points = np.column_stack( (grid_x.ravel(), grid_y.ravel()) )
    # grid_values = np.full(grid_points.shape[0], np.nan, dtype=np.float32)
    
    # Convert to Dask arrays if not already
    # obs_points_da = da.from_array(obs_points, chunks=(chunk_size, 2))
    # obs_values_da = da.from_array(obs_values, chunks=(chunk_size,))
    grid_points_da = da.from_array(grid_points, chunks=(chunk_size, 2))
    
    radius2 = radius**2
    
    def interpolate_point(gp):
        """Interpolate a single grid point (NumPy inside Dask map_blocks)."""
        dx = obs_points[:, 0] - gp[0] #should this be obs_points_da?
        dy = obs_points[:, 1] - gp[1]
        dist2 = dx**2 + dy**2

        mask = dist2 <= radius2
        if not np.any(mask):
            return np.nan

        weights = (radius2 - dist2[mask]) / (radius2 + dist2[mask])
        weights = np.clip(weights, 0, None)

        if np.sum(weights) > 0:
            return np.sum(weights * obs_values[mask]) / np.sum(weights)
        else:
            return np.nan
        
    # Apply interpolation to each grid point lazily
    interpolated_da = grid_points_da.map_blocks(
        lambda block: np.array([interpolate_point(gp) for gp in block]),
        dtype=float,
        drop_axis=1
    )

    # Reshape to grid
    interpolated_da = interpolated_da.reshape(grid_x.shape)
    
    # Compute with progress bar
    with ProgressBar():
        interpolated = interpolated_da.compute()

    return interpolated




def longest_continuous_branch(multilines):
    if not isinstance(multilines, MultiLineString):
        raise TypeError("Input must be a shapely MultiLineString")

    # Build a graph where nodes are coordinates and edges are line segments
    G = nx.Graph()
    for line in multilines.geoms:
        coords = list(line.coords)
        for i in range(len(coords) - 1):
            p1, p2 = coords[i], coords[i + 1]
            length = LineString([p1, p2]).length
            G.add_edge(p1, p2, weight=length)

    # Find connected components
    longest_path = []
    max_length = 0.0

    for component in nx.connected_components(G):
        subgraph = G.subgraph(component)
        # Find all pairs shortest paths (weighted by length)
        for u in subgraph.nodes:
            lengths, paths = nx.single_source_dijkstra(subgraph, u, weight="weight")
            for v, length in lengths.items():
                if length > max_length:
                    max_length = length
                    longest_path = paths[v]

    # Convert longest path to a LineString
    if longest_path:
        return LineString(longest_path), max_length
    else:
        return None, 0.0


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
    
    # # this does not work and idk why
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





#%%


def reduce_floats(variables, nbits=32):
    # variables: list of floats and/or floating point arrays
    # nbits: scalar, default 32
    
    reduced_variables = []
    
    for v in variables:
        if nbits == 32:
            reduced_variables.append( v.astype(np.float32) )
        elif nbits == 16:
            reduced_variables.append( v.astype(np.float16) )
    
    return reduced_variables


































