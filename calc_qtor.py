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
# from glob import glob
from scipy.ndimage import gaussian_filter
import os
from os.path import exists
from matplotlib.ticker import MultipleLocator


### STEPS
# 1. QLCS object identification
# 2. QLCS leading line identification
# 3. Storm motion estimation
# 4. Leading line shear analysis
# 5. Local tortuosity
# 6. Calculate QTor


# NEXRAD fields: 'differential_reflectivity', 'reflectivity', 'cross_correlation_ratio', 'clutter_filter_power_removed', 'spectrum_width', 'differential_phase', 'velocity'

# ECCC radars: geoJSON files?
# From IWPro source code
# https://cssl.inwx.ca/static/js/main.2f8e1c6f.chunk.js
# 'H5': {
#        'tilt': ['1', '2', '3', '4'],
#        'ref_tilt': ['1', '2', '3', '4'],
#        'products': ['Reflectivity', 'Velocity', 'Storm\x20Relative\x20Velocity', 'Correlation\x20Coefficient', 'Hail\x20Differential\x20Reflectivity', 'Tornado\x20Debris\x20Signature', 'Differential\x20Reflectivity', 'Specific\x20Differential\x20Phase', 'Spectrum\x20Width'],
#        'ref_products': ['ref', 'vel', 'srv', 'cc', 'hdr', 'tds', 'dr', 'sdp', 'sw']
#        },
# 'H5_CANADIAN': {
#                 'tilt': ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17'],
#                 'ref_tilt': ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17'],
#                 'products': ['Reflectivity', 'Velocity', 'Storm\x20Relative\x20Velocity', 'Correlation\x20Coefficient', 'Hail\x20Differential\x20Reflectivity', 'Tornado\x20Debris\x20Signature', 'Differential\x20Reflectivity', 'Specific\x20Differential\x20Phase', 'Spectrum\x20Width'],
#                 'ref_products': ['ref', 'vel', 'srv', 'cc', 'hdr', 'tds', 'dr', 'sdp', 'sw']



# 'https://'['concat']('instantweather-reports', '.s3.amazonaws.com/')['concat'](String(_0x3d991f || '')['split']('/')['map'](function(_0x3b444b)
# 'S3_CA_PATH': '/alerts/ca',
# 'S3_US_PATH': '/alerts/us-v3',

#%% QLCS object criteria
# Britt et al. 2026 https://doi.org/10.1175/WAF-D-25-0092.1

# Britt et al. 2024 QLCS criteria https://doi.org/10.1175/WAF-D-23-0106.1
min_cref = 40 #composite dbz threshold to get initial objects
max_ref = 45 #max object dbz threshold
min_area = 54 #storm object area threshold [km^2]
merge_dist = 12 #proximity threshold [km]
min_length_1 = 100 #length threshold [km]
min_length_2 = 150 #length threshold [km]
min_ecc_1 = 0.85 #eccentricity threshold for length>100 km
min_ecc_2 = 0.74 #eccentricity threshold for length>150 km


# Potvin et al. 2022 QLCS criteria https://doi.org/10.1175/JTECH-D-21-0141.1



# Smith et al. 2012 QLCS criteria https://doi.org/10.1175/WAF-D-11-00115.1
min_cref = 35 #dbz threshold
min_length = 100 #length threshold for 35+ dBZ area [km]
min_aspect = 3 #aspect ratio threshold







#%% 


























