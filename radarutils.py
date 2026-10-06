# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 12:10:53 2026

@author: mschne28
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
import cfgrib
import pandas as pd
import h5py as h5
# import sklearn
from glob import glob
import os
from os.path import exists
from matplotlib.ticker import MultipleLocator
from datetime import datetime

from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from skimage import measure,morphology,filters
# from metpy.interpolate import interpolate_to_points
from scipy.ndimage import gaussian_filter1d
import math


# NEXRAD fields: 'differential_reflectivity', 'reflectivity', 'cross_correlation_ratio', 'clutter_filter_power_removed', 'spectrum_width', 'differential_phase', 'velocity'
# ECCC fields: 'reflectivity_horizontal', 'total_power_horizontal', 'cross_correlation_ratio', 'differential_phase', 'differential_reflectivity', 'specific_differential_phase', 'velocity_horizontal'



#%% ECCC radar functions

# Read ECCC radar xata - tbd (to be developed) once I have the data
def read_eccc_hdf5(filename, grid_x, grid_y, max_rng=300, radius=3):
    # Text header block needs t be removed from HDF5 files (first 34 bytes)
    fp = os.path.dirname(filename)
    file_ext_ind = filename.find('.h')
    filename_stripped = filename[:file_ext_ind] + "_stripped" + filename[file_ext_ind:]
    h5_files = [f.replace("\\", "/") for f in glob(fp + "/*")]
    
    if filename_stripped not in h5_files:
        header_length = 34
        chunk_size = 1024*1024
        with open(filename, 'rb') as src, open(filename_stripped, 'wb') as dst:
            src.seek(header_length)
            while True:
                chunk = src.read(chunk_size)
                if not chunk:
                    break
                dst.write(chunk)
    
    # Read file with header removed
    radar = pyart.aux_io.read_odim_h5(filename_stripped)
    
    radar.init_gate_x_y_z()
    radar.init_gate_longitude_latitude()
    radar.init_gate_altitude()
    radar.init_rays_per_sweep()
    
    return radar




'''
### Resource: https://data.eol.ucar.edu/zinc/file/download/54E3335CD84EE/OPERA2014_O4_ODIM_H5-v2.2.pdf ###

# Entire HDF structure (f, level 0) --> subgroups 'dataset1'-'dataset17' + metadata 'where'], 'what', and 'how'
f['where'] --> radar location
        attributes ['height', 'lat', 'lon']
f['what'] --> 
        attributes ['date', 'object', 'source', 'time', 'version']
f['how'] --> radar operating parameters
        attributes ['L_UPDATED_TASK', 'RXlossH', 'RXlossV', 'TXlossH', 'TXlossV', 'TXtype', '_creator_program', '_orig_file_format', '_orig_sensor_id', '_orig_sensor_name',
                    'antgainH', 'antgainV', 'beamwH', 'beamwV', 'beamwidth', 'frequency', 'polmode', 'poltype', 'radomelossH', 'radomelossV', 'scan_count', 'simulated',
                    'software', 'sw_version', 'system', 'task', 'time_accuracy_downgrade']

# Groups (level 1) 'dataset1'-'dataset17' --> these are the individual sweeps
# Note: ECCC radar volumes start at the highest tilt, not the lowest!
 ['datasetn'] --> subgroups 'data1'-'data10' + metadata groups 'what', 'where', and 'how'
 ['datasetn']['where']
       attributes ['a1gate', 'elangle', 'nbins', 'nrays', 'rscale', 'rstart'] --> access via f['datasetn']['where'].attrs[key]
 ['datasetn']['what']
       attributes ['enddate', 'endtime', 'product', 'startdate', 'starttime']
 ['datasetn']['how']
       attributes ['CSR', 'LOG', 'NEZH', 'NEZH_A', 'NEZV', 'NEZV_A', 'NI', 'SQI', 'TXcalpowkwH', 'TXcalpowkwV', 'TXpower', 'Vsamples', 'anglesync', 'anglesyncRes',
                    'antspeed', 'astart', 'avgpwr', 'azangles', 'azmethod', 'base_1km_hc', 'base_1km_vc', 'binmethod', 'binmethod_avg', 'clutterType', 'comment',
                    'dataflag', 'dielectic_factor', 'dual-pol_TXpower', 'dual-pol_avgpwr', 'dual-pol_peakpwr', 'elangles', 'highprf', 'lowprf', 'malfunc', 'numpulses',
                    'peakpwr', 'phasediff', 'pol_of_txpower', 'pulsewidth', 'radar_msg', 'radconstH', 'radconstV', 'scan_index', 'single-pol_TXpower', 'single-pol_avgpwr',
                    'single-pol_peakpwr', 'startT', 'startazA', 'startelA', 'stopazA', 'stopelA', 'zcalH', 'zcalV', 'zdrcal']

# Sub-groups (level 2) 'data1'-'data10' --> these are the different variables
['datasetn']['datan'] --> subgroups 'data', 'what'
['datasetn']['datan']['data'][:] --> actual data array
['datasetn']['datan']['what'].attrs = ['gain', 'nodata', 'offset', 'quantity', 'undetect']
    Variable names stored in ['datasetn']['datan']['what'].attrs['quantity']
    data1  = DBZH    - Corrected horizontal reflectivity factor (dBZ)
    data2  = RHOHV   - Correlation coefficient (0-1)
    data3  = UPHIDP  - Uncorrected PHIDP
    data4  = WRADH   - Spectrum width (horizontal)
    data5  = PHIDP   - Differential phase (deg)
    data6  = ZDR     - Differential reflectivity (dBZ? or dB?)
    data7  = KDP     - Specific differential phase (deg/km)
    data8  = SQIH    - Horizontal signal quality index (0-1)
    data9  = VRADH   - Velocity (horizontal) (m/s)
    data10 = TH      - Total uncorrected horizontal reflectivity factor (dBZ)


volume_date = f['what'].attrs['date'].decode('utf-8')
volume_time = f['what'].attrs['time'].decode('utf-8')
time = datetime.strptime(volume_date + volume_time, '%Y%m%d%H%M%S')

radar_height = f['where'].attrs['height']
radar_lat = f['where'].attrs['lat']
radar_lon = f['where'].attrs['lon']

for i in range(17):
    dat = f[f"dataset{i+1}"]
    
    elev_angle = dat['where'].attrs['elangle']
    nbins = dat['where'].attrs['nbins']
    nrays = dat['where'].attrs['nrays']
    dr = dat['where'].attrs['rscale'] #range resolution (m)
    
    az_angles = dat['how'].attrs['azangles']
    
    dbz = dat['data1']['data'][:]
    dbz_gain = dat['data1']['what'].attrs['gain'] # a in y=ax+b used to convert to unit?
    dbz_offset = dat['data1']['what'].attrs['offset'] # b in y=ax+b used to convert to unit?
    dbz_name = dat['data1']['what'].attrs['quantity'].decode('utf-8')
    
    vel = dat['data9']['data'][:]
    zdr = dat['data6']['data'][:]
    rhohv = dat['data2']['data'][:]
    kdp = dat['data7']['data'][:]
    sw = dat['data4']['data'][:]


'''


# def read_eccc_geojson(filename, radar_name, max_rng=300):
#     df = pd.read_csv('C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/ECCC_radar_locations.csv',
#                      sep=",", header=0, usecols=["Call sign", "Latitude", "Longitude"], index_col="Call sign")
#     radar_lat = df.loc[radar_name]['Latitude']
#     radar_lon = df.loc[radar_name]['Longitude']

#     return radar_lat, radar_lon




# Read NEXRAD data
def read_nexrad(filename, max_rng=300):
    radar = pyart.io.read(filename)
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_below('cross_correlation_ratio', 0.9)
    
    kdp, phidp = pyart.retrieve.kd_maesaka(radar, gatefilter=gatefilter, phidp_field='differential_phase')
    radar.add_field('specific_differential_phase', kdp)
    
    return radar




