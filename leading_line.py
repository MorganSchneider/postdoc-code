# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 16:28:31 2026

@author: mschne28
"""

import matplotlib.pyplot as plt
import numpy as np
import shapely
from shapely import Polygon
import pygeoops as pg
from skimage import measure,morphology,filters
import pickle


#%% 


import shapely
from shapely.geometry import Polygon,LineString,MultiLineString
from shapely.plotting import plot_polygon
import pygeoops as pg


fp = 'C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/qlcs_tornado_outbreaks/'

dbfile = open(fp+'qlcs_object_test.pkl', 'rb')
obj = pickle.load(dbfile)
qlcs_obj = obj['qlcs_obj']
cref_smooth = obj['cref_smooth']
xm = obj['xm']
ym = obj['ym']
x_centroid = obj['x_centroid']
y_centroid = obj['y_centroid']
dbfile.close()

dbfile = open(fp+'leading_line_test.pkl', 'rb')
ll = pickle.load(dbfile)
dbfile.close()



qlcs_obj_smooth = filters.gaussian(qlcs_obj, sigma=5/3)

thres = np.percentile(qlcs_obj_smooth[(qlcs_obj_smooth>0)], 60) #could just use 5e-20? idk if the limits are universal
# thres = np.max(qlcs_obj_smooth) / 2
# thres = 5e-20
obj_smooth_bin = np.where(qlcs_obj_smooth>thres, 1, 0)

contours = measure.find_contours(obj_smooth_bin.astype(bool), level=0.99, fully_connected='low')
contours_int = [np.round(contours[n]).astype(int) for n in range(len(contours))]

# contour = np.concatenate(contours)
ind = np.asarray([len(contours[n]) for n in range(len(contours))])
contour = contours[np.argmax(ind)]
obj_contour_int = np.round(contour).astype(int)
obj_contour = np.array([ [xm[j,i], ym[j,i]] for j,i in zip(obj_contour_int[:,0], obj_contour_int[:,1]) ])


polygon = Polygon(obj_contour)
# geom = shapely.from_wkt(polygon.wkt)
# shapely.remove_repeated_points(geom)


fig,ax = plt.subplots(1, 1, figsize=(8,6))
c = ax.pcolormesh(xm, ym, cref_smooth, vmin=0, vmax=70, cmap='gist_ncar')
cb = plt.colorbar(c, ax=ax, extend='both')
plot_polygon(polygon, ax=ax, add_points=False, facecolor=None, edgecolor='k', linewidth=1)
ax.scatter(x_centroid, y_centroid, marker='.', s=30, c='k')
ax.set_title(f"QLCS object smoothed, binarized")
plt.show()
#%%

cl = pg.centerline(polygon, min_branch_length=50)


cl




fig,ax = plt.subplots(1, 1, figsize=(8,6))
c = ax.pcolormesh(xm, ym, cref_smooth, vmin=0, vmax=70, cmap='gist_ncar')
cb = plt.colorbar(c, ax=ax, extend='both')
plot_polygon(polygon, ax=ax, add_points=False, facecolor=None, edgecolor='k', linewidth=1)

ax.scatter(x_centroid, y_centroid, marker='.', s=30, c='k')
ax.set_title(f"QLCS object smoothed, binarized")
plt.show()







