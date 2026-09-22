# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 11:15:55 2026

@author: mschne28
"""

# Calculate QTor


from QTORutils import *


### STEPS
# 1. QLCS object identification DONEEEEEE
# 2. QLCS leading line identification
# 3. Storm motion estimation
# 4. Leading line shear analysis
# 5. Local tortuosity
# 6. Calculate QTor


# NEXRAD fields: 'differential_reflectivity', 'reflectivity', 'cross_correlation_ratio', 'clutter_filter_power_removed', 'spectrum_width', 'differential_phase', 'velocity'


### Functions ###
# cressman_interpolation_dask
# find_qlcs_objects
# distance
# get_storm_motion
# get_centroid
# get_leading_line
# get_inflow_polygon
# get_shear
# InterpolateToHeightAboveGround
# 


### Functions to write
# get_line_normal_shear
# get_line_parallel_shear
# get_local_tortuosity
# calc_QTor





#%% Open file

fp = 'C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/qlcs_tornado_outbreaks/'

radar_name = 'KBUF'
yyyyt = 2026
mmt = 8
ddt = 2
timestr = '162824'
# timestr = '161342'
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
gate_lat = cradar.extract_sweeps([0]).gate_latitude['data']
gate_lon = cradar.extract_sweeps([0]).gate_longitude['data']

#%% Cressman filtering



max_rng = 300 #set max range, km
hres = 3 #gridded resolution, km
radius = 3

# pts = np.array([gate_x[:,(rng<=max_rng)].ravel(), gate_y[:,(rng<=max_rng)].ravel()]).transpose()
# vals = cref[:,(rng<=max_rng)].ravel().data
# xm,ym = np.meshgrid(np.arange(-max_rng, max_rng+hres, hres), np.arange(-max_rng, max_rng+hres, hres))
# xi = np.array([xm.ravel(), ym.ravel()]).transpose()
# cref_interp_flat = interpolate_to_points(pts, vals, xi, interp_type='linear', search_radius=hres)
# cref_interp = cref_interp_flat.reshape(xm.shape)
# cref_smooth = filters.gaussian(cref_interp, sigma=0.8)


gx = gate_x[:,(rng<=max_rng)].astype(np.float32)
gy = gate_y[:,(rng<=max_rng)].astype(np.float32)
vals = cref[:,(rng<=max_rng)].data
xm,ym = np.meshgrid(np.arange(-max_rng, max_rng+hres, hres), np.arange(-max_rng, max_rng+hres, radius))
cref_smooth = cressman_interpolation_dask(gx, gy, vals, xm, ym, radius, chunk_size=5000)


# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
# ax.set_title('Cressman interpolated cref')
# plt.show()


glat = gate_lat[:,(rng<=max_rng)].astype(np.float32)
glon = gate_lon[:,(rng<=max_rng)].astype(np.float32)
lonm,latm = np.meshgrid(np.linspace(np.min(glon), np.max(glon), len(xm)), np.linspace(np.min(glat), np.max(glat), len(xm)))

#%% Retrieve QLCS object IDs


qlcs_labels, qlcs_regions = get_qlcs_objects(cref_smooth, hres)

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

# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, qlcs_obj, 'dbz', ax, datalims=[0,1], cmap='HomeyerRainbow')
# ax.set_title('Final QLCS object\n Object with highest dBZ ')
# plt.show()



#%% Storm motion estimation - get data from previous and next volumes


filenames = glob(fp+f"{radar_name}/{radar_name}{yyyyt}{mmt:02.0f}{ddt:02.0f}_*_V06.ar2v")
ind = [i for i in range(len(filenames)) if timestr in filenames[i]][0]


# Previous volume
radar1 = pyart.io.read(filenames[ind-1])
gatefilter = pyart.filters.GateFilter(radar1)
gatefilter.exclude_transition()
gatefilter.exclude_below('cross_correlation_ratio', 0.9)
cradar_prev = pyart.retrieve.composite_reflectivity(radar1, field='reflectivity', gatefilter=gatefilter)
cref_prev = cradar_prev.fields['composite_reflectivity']['data']

gate_x = cradar_prev.extract_sweeps([0]).gate_x['data']/1000
gate_y = cradar_prev.extract_sweeps([0]).gate_y['data']/1000

print('Previous volume')
gx = gate_x[:,(rng<=max_rng)].astype(np.float32)
gy = gate_y[:,(rng<=max_rng)].astype(np.float32)
vals = cref_prev[:,(rng<=max_rng)].data
cref_smooth_prev = cressman_interpolation_dask(gx, gy, vals, xm, ym, radius, chunk_size=5000)




# Next volume
radar2 = pyart.io.read(filenames[ind+1])
gatefilter = pyart.filters.GateFilter(radar2)
gatefilter.exclude_transition()
gatefilter.exclude_below('cross_correlation_ratio', 0.9)
cradar_next = pyart.retrieve.composite_reflectivity(radar2, field='reflectivity', gatefilter=gatefilter)
cref_next = cradar_next.fields['composite_reflectivity']['data']

gate_x = cradar_next.extract_sweeps([0]).gate_x['data']/1000
gate_y = cradar_next.extract_sweeps([0]).gate_y['data']/1000

print('Next volume')
gx = gate_x[:,(rng<=max_rng)].astype(np.float32)
gy = gate_y[:,(rng<=max_rng)].astype(np.float32)
vals = cref_next[:,(rng<=max_rng)].data
cref_smooth_next = cressman_interpolation_dask(gx, gy, vals, xm, ym, radius, chunk_size=5000)



# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, cref_smooth_prev, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
# ax.set_title('Cressman interpolated cref, previous volume')
# plt.show()

# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, cref_smooth_next, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
# ax.set_title('Cressman interpolated cref, next volume')
# plt.show()



#%% Find QLCS objects in previous and next volumes

qlcs_labels_prev, qlcs_regions_prev = get_qlcs_objects(cref_smooth_prev, hres)
qlcs_labels_next, qlcs_regions_next = get_qlcs_objects(cref_smooth_next, hres)



#%% Get centroids from analysis scan + previous and next scans

# jb1,ib1,jb2,ib2 = qlcs_region.bbox
# bbox = [[xm[jb1,ib1], ym[jb1,ib1]], [xm[jb1,ib2-1], ym[jb1,ib2-1]], [xm[jb2-1,ib1], ym[jb2-1,ib1]], [xm[jb2-1,ib2-1], ym[jb2-1,ib2-1]]]

j_centroid,i_centroid = qlcs_region.centroid
jc1,ic1,jc2,ic2 = np.floor(j_centroid), np.floor(i_centroid), np.ceil(j_centroid), np.ceil(i_centroid)
j1,i1,j2,i2 = int(jc1), int(ic1), int(jc2), int(ic2)
x_centroid = (1 - abs(i_centroid-ic1))*xm[j1,i1] + (1 - abs(i_centroid-ic2))*xm[j1,i2]
y_centroid = (1 - abs(j_centroid-jc1))*ym[j1,i1] + (1 - abs(j_centroid-jc2))*ym[j2,i1]
centroid = [x_centroid, y_centroid]


#% Get centroids from previous and next scan

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
        j_centroid,i_centroid = qlcs_regions_prev[n].centroid
        jc1,ic1,jc2,ic2 = np.floor(j_centroid), np.floor(i_centroid), np.ceil(j_centroid), np.ceil(i_centroid)
        j1,i1,j2,i2 = int(jc1), int(ic1), int(jc2), int(ic2)
        xc[n] = (1 - abs(i_centroid-ic1))*xm[j1,i1] + (1 - abs(i_centroid-ic2))*xm[j1,i2]
        yc[n] = (1 - abs(j_centroid-jc1))*ym[j1,i1] + (1 - abs(j_centroid-jc2))*ym[j2,i1]
        
        centroid_dist2[n] = (xc[n] - x_centroid)**2 + (yc[n] - y_centroid)**2
        
    n_closest = np.argmin(centroid_dist2)
    
    qlcs_obj_prev[(qlcs_labels_prev == n_closest+1)] = 1
    # qlcs_obj_prev[(qlcs_labels_prev == regs[n_closest])] = 1
    qlcs_region_prev = qlcs_regions_prev[n_closest]
    x_centroid_prev = xc[n_closest]
    y_centroid_prev = yc[n_closest]
else:
    j_centroid,i_centroid = qlcs_regions_prev[0].centroid
    jc1,ic1,jc2,ic2 = np.floor(j_centroid), np.floor(i_centroid), np.ceil(j_centroid), np.ceil(i_centroid)
    j1,i1,j2,i2 = int(jc1), int(ic1), int(jc2), int(ic2)
    x_centroid_prev = (1 - abs(i_centroid-ic1))*xm[j1,i1] + (1 - abs(i_centroid-ic2))*xm[j1,i2]
    y_centroid_prev = (1 - abs(j_centroid-jc1))*ym[j1,i1] + (1 - abs(j_centroid-jc2))*ym[j2,i1]
    
    qlcs_obj_prev = qlcs_labels_prev
    


if len(qlcs_regions_next) > 1:
    xc = np.zeros((len(qlcs_regions_next),))
    yc = np.zeros((len(qlcs_regions_next),))
    centroid_dist2 = np.zeros((len(qlcs_regions_next),))
    regs = np.unique(qlcs_labels_next[(qlcs_labels_next>0)])
    for n in range(len(centroid_dist2)):
        j_centroid,i_centroid = qlcs_regions_next[n].centroid
        jc1,ic1,jc2,ic2 = np.floor(j_centroid), np.floor(i_centroid), np.ceil(j_centroid), np.ceil(i_centroid)
        j1,i1,j2,i2 = int(jc1), int(ic1), int(jc2), int(ic2)
        xc[n] = (1 - abs(i_centroid-ic1))*xm[j1,i1] + (1 - abs(i_centroid-ic2))*xm[j1,i2]
        yc[n] = (1 - abs(j_centroid-jc1))*ym[j1,i1] + (1 - abs(j_centroid-jc2))*ym[j2,i1]
        
        centroid_dist2[n] = (xc[n] - x_centroid)**2 + (yc[n] - y_centroid)**2
        
    n_closest = np.argmin(centroid_dist2)
    
    qlcs_obj_next[(qlcs_labels_next == n_closest+1)] = 1
    # qlcs_obj_next[(qlcs_labels_next == regs[n_closest])] = 1
    qlcs_region_next = qlcs_regions_next[n_closest]
    x_centroid_next = xc[n_closest]
    y_centroid_next = yc[n_closest]
else:
    j_centroid,i_centroid = qlcs_regions_next[0].centroid
    jc1,ic1,jc2,ic2 = np.floor(j_centroid), np.floor(i_centroid), np.ceil(j_centroid), np.ceil(i_centroid)
    j1,i1,j2,i2 = int(jc1), int(ic1), int(jc2), int(ic2)
    x_centroid_next = (1 - abs(i_centroid-ic1))*xm[j1,i1] + (1 - abs(i_centroid-ic2))*xm[j1,i2]
    y_centroid_next = (1 - abs(j_centroid-jc1))*ym[j1,i1] + (1 - abs(j_centroid-jc2))*ym[j2,i1]
    
    qlcs_obj_next = qlcs_labels_next


#%% Estimate storm motion from centroid displacement


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



#%% Initial leading line identification using storm motion

# E-W search
x_ll_u = []
y_ll_u = []
idx_ll_u = []
jdy_ll_u = []
lat_ll_u = []
lon_ll_u = []

for jdy in range(xm.shape[0]):
    if np.any(qlcs_obj[jdy,:] > 0):
        if u_sm > 0:
            idx = np.where(qlcs_obj[jdy,:]>0)[0][-1]
        elif u_sm < 0:
            idx = np.where(qlcs_obj[jdy,:]>0)[0][0]
        
        x_ll_u.append(xm[jdy,idx])
        y_ll_u.append(ym[jdy,idx])
        idx_ll_u.append(idx)
        jdy_ll_u.append(jdy)
        lat_ll_u.append(latm[jdy,idx])
        lon_ll_u.append(lonm[jdy,idx])
        

# N-S search
x_ll_v = []
y_ll_v = []
idx_ll_v = []
jdy_ll_v = []
lat_ll_v = []
lon_ll_v = []

for idx in range(ym.shape[1]):
    if np.any(qlcs_obj[:,idx] > 0):
        if v_sm > 0:
            jdy = np.where(qlcs_obj[:,idx]>0)[0][-1]
        elif v_sm < 0:
            jdy = np.where(qlcs_obj[:,idx]>0)[0][0]
        
        x_ll_v.append(xm[jdy,idx])
        y_ll_v.append(ym[jdy,idx])
        idx_ll_v.append(idx)
        jdy_ll_v.append(jdy)
        lat_ll_v.append(latm[jdy,idx])
        lon_ll_v.append(lonm[jdy,idx])


#%% Save leading line to pickles

ll_zonal = {'x':x_ll_u, 'y':y_ll_u, 'i':idx_ll_u, 'j':jdy_ll_u, 'lat':lat_ll_u, 'lon':lon_ll_u}
ll_meridional = {'x':x_ll_v, 'y':y_ll_v, 'i':idx_ll_v, 'j':jdy_ll_v, 'lat':lat_ll_v, 'lon':lon_ll_v}
data = {'zonal':ll_zonal, 'meridional':ll_meridional, 'u_sm':u_sm, 'v_sm':v_sm}

dbfile = open(fp+'leading_line_test.pkl', 'wb')
pickle.dump(data, dbfile)
dbfile.close()



data = {'xm':xm, 'ym':ym, 'latm':latm, 'lonm':lonm, 'cref_smooth':cref_smooth, 'qlcs_obj':qlcs_obj, 'x_centroid':x_centroid, 'y_centroid':y_centroid, 'qlcs_region':qlcs_region}
dbfile = open(fp+'qlcs_object_test.pkl', 'wb')
pickle.dump(data, dbfile)
dbfile.close()




#%% Find centerline


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



# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, qlcs_obj_smooth, 'dbz', ax, datalims=[0,np.max(qlcs_obj_smooth)], cmap='HomeyerRainbow')
# ax.scatter(x_centroid, y_centroid, marker='.', s=30, c='k')
# ax.set_title(f"QLCS object smoothed")
# plt.show()

# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, obj_smooth_bin, 'dbz', ax, datalims=[0,1], cmap='HomeyerRainbow')
# ax.scatter(x_centroid, y_centroid, marker='.', s=30, c='k')
# for j,i in zip(obj_contour_int[:,0], obj_contour_int[:,1]):
#     ax.scatter(xm[j,i], ym[j,i], marker='.', s=1, c='yellow')
# ax.set_title(f"QLCS object smoothed, binarized")
# plt.show()




polygon = Polygon(obj_contour)

centerline = get_centerline(polygon)
trimmed_centerline = substring(centerline, 0.05*centerline.length, 0.95*centerline.length)

# centerline = Centerline(polygon)
# if True:
#     new_centerline, centerline_length = longest_continuous_branch(centerline.geometry)
#     centerline = new_centerline
#     trimmed_centerline = substring(new_centerline, 0.05*centerline_length, 0.95*centerline_length)

x_cl,y_cl = trimmed_centerline.xy

#%

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
plot_polygon(polygon, ax=ax, add_points=False, facecolor=None, edgecolor='k', linewidth=1)
ax.plot(x_cl, y_cl, 'k', linewidth=1.5)
ax.scatter(x_centroid, y_centroid, marker='.', s=30, c='k')
ax.set_title(f"QLCS polygon with centerline")
plt.show()



if True:
    data_cl = {'x':x_cl, 'y':y_cl, 'centerline':centerline, 'trimmed_centerline':trimmed_centerline}
    dbfile = open(fp+'centerline_test.pkl', 'wb')
    pickle.dump(data_cl, dbfile)
    dbfile.close()



#%% QC and smooth leading line points
#   Find closest centerline point to each leading line point,
#   create vector between centerline and leading line points,
#   exclude any leading line points misaligned with expected SM



dbfile = open(fp+'leading_line_test.pkl', 'rb')
ll = pickle.load(dbfile)
u_sm = ll['u_sm']
v_sm = ll['v_sm']
x_zonal = ll['zonal']['x']
y_zonal = ll['zonal']['y']
x_merid = ll['meridional']['x']
y_merid = ll['meridional']['y']
lat_zonal = ll['zonal']['lat']
lon_zonal = ll['zonal']['lon']
lat_merid = ll['meridional']['lat']
lon_merid = ll['meridional']['lon']
dbfile.close()


dbfile = open(fp+'centerline_test.pkl', 'rb')
cl = pickle.load(dbfile)
x_cl = cl['x']
y_cl = cl['y']
dbfile.close()


cl_points = np.column_stack( (x_cl, y_cl))
zonal_points = np.column_stack( (x_zonal, y_zonal))
merid_points = np.column_stack( (x_merid, y_merid))
zonal_lonlat = np.column_stack( (lon_zonal, lat_zonal))
merid_lonlat = np.column_stack( (lon_merid, lat_merid))

tree_cl = KDTree(cl_points)


# Zonal
dist_zonal,inds_zonal = tree_cl.query(zonal_points) # Find nearest centerline point in tree_cl to each zonal point

# Meridional
dist_merid,inds_merid = tree_cl.query(merid_points) # Find nearest centerline point in tree_cl to each meridional point

cl_points_matched_zonal = np.array([cl_points[inds_zonal[i]] for i in range(len(inds_zonal))])
cl_points_matched_merid = np.array([cl_points[inds_merid[i]] for i in range(len(inds_merid))])

vector_zonal = zonal_points - cl_points_matched_zonal
vector_merid = merid_points - cl_points_matched_merid






# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
# plot_polygon(polygon, ax=ax, add_points=False, facecolor=None, edgecolor='k', linewidth=1)
# ax.plot(x_cl, y_cl, 'k', linewidth=1.5)
# ax.scatter(x_zonal, y_zonal, marker='.', s=1, c='red')
# ax.plot(np.column_stack((cl_points_matched_zonal[:,0], x_zonal)).transpose(), np.column_stack((cl_points_matched_zonal[:,1], y_zonal)).transpose(), '-r', linewidth=0.75)
# ax.set_title(f"QLCS polygon with centerline and zonal points")
# ax.set_xlim([-200,100])
# ax.set_ylim([-100,200])
# plt.show()




# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
# plot_polygon(polygon, ax=ax, add_points=False, facecolor=None, edgecolor='k', linewidth=1)
# ax.plot(x_cl, y_cl, 'k', linewidth=1.5)
# ax.scatter(x_merid, y_merid, marker='.', s=1, c='red')
# ax.plot(np.column_stack((cl_points_matched_merid[:,0], x_merid)).transpose(), np.column_stack((cl_points_matched_merid[:,1], y_merid)).transpose(), '-r', linewidth=0.75)
# ax.set_title(f"QLCS polygon with centerline and meridional points")
# ax.set_xlim([-200,100])
# ax.set_ylim([-100,200])
# plt.show()



if abs(u_sm) > abs(v_sm):
    # zonal_points_filtered = zonal_points[((vector_zonal[:,0]/u_sm > 0) & (abs(vector_zonal[:,0])>abs(vector_zonal[:,1]))),:]
    # merid_points_filtered = merid_points[((vector_merid[:,0]/u_sm > 0) & (abs(vector_merid[:,0])>abs(vector_merid[:,1]))),:]
    # cl_points_filtered_zonal = cl_points_matched_zonal[((vector_zonal[:,0]/u_sm > 0) & (abs(vector_zonal[:,0])>abs(vector_zonal[:,1]))),:]
    zonal_points_filtered = zonal_points[(vector_zonal[:,0]/u_sm > 0),:]
    merid_points_filtered = merid_points[(vector_merid[:,0]/u_sm > 0),:]
    zonal_lonlat_filtered = zonal_lonlat[(vector_zonal[:,0]/u_sm > 0),:]
    merid_lonlat_filtered = merid_lonlat[(vector_merid[:,0]/u_sm > 0),:]
    idx_zonal_filtered = np.asarray(ll['zonal']['i'])[(vector_zonal[:,0]/u_sm > 0)]
    jdy_zonal_filtered = np.asarray(ll['zonal']['j'])[(vector_zonal[:,0]/u_sm > 0)]
    idx_merid_filtered = np.asarray(ll['meridional']['i'])[(vector_merid[:,0]/u_sm > 0)]
    jdy_merid_filtered = np.asarray(ll['meridional']['j'])[(vector_merid[:,0]/u_sm > 0)]
elif abs(v_sm) > abs(u_sm):
    # zonal_points_filtered = zonal_points[((vector_zonal[:,1]/v_sm > 0) & (abs(vector_zonal[:,1])>abs(vector_zonal[:,0]))),:]
    # merid_points_filtered = merid_points[((vector_merid[:,1]/v_sm > 0) & (abs(vector_merid[:,1])>abs(vector_merid[:,0]))),:]
    zonal_points_filtered = zonal_points[(vector_zonal[:,1]/v_sm > 0),:]
    merid_points_filtered = merid_points[(vector_merid[:,1]/v_sm > 0),:]
    zonal_lonlat_filtered = zonal_lonlat[(vector_zonal[:,1]/v_sm > 0),:]
    merid_lonlat_filtered = merid_lonlat[(vector_merid[:,1]/v_sm > 0),:]
    idx_zonal_filtered = np.asarray(ll['zonal']['i'])[(vector_zonal[:,1]/v_sm > 0)]
    jdy_zonal_filtered = np.asarray(ll['zonal']['j'])[(vector_zonal[:,1]/v_sm > 0)]
    idx_merid_filtered = np.asarray(ll['meridional']['i'])[(vector_merid[:,1]/v_sm > 0)]
    jdy_merid_filtered = np.asarray(ll['meridional']['j'])[(vector_merid[:,1]/v_sm > 0)]





# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
# plot_polygon(polygon, ax=ax, add_points=False, facecolor=None, edgecolor='k', linewidth=1)
# ax.plot(x_cl, y_cl, 'k', linewidth=1.5)
# ax.scatter(zonal_points_filtered[:,0], zonal_points_filtered[:,1], marker='.', s=2, c='r')
# ax.scatter(merid_points_filtered[:,0], merid_points_filtered[:,1], marker='.', s=2, c='b')
# # ax.plot(np.column_stack((cl_points_matched_zonal[:,0], x_zonal)).transpose(), np.column_stack((cl_points_matched_zonal[:,1], y_zonal)).transpose(), '-r', linewidth=0.75)
# ax.set_title(f"QLCS polygon with centerline and zonal points")
# ax.set_xlim([-200,100])
# ax.set_ylim([-100,200])
# plt.show()






def distance(point1, point2):
    return math.hypot(point2[0] - point1[0], point2[1] - point1[1])




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
query_point = [-300, 300]
points_tmp = [item for item in ll_points_cat]
lonlat_tmp = [item for item in ll_lonlat_cat]
is_sorted = [False] * len(ll_points_cat)

tree_ll = KDTree(np.asarray(ll_points_cat))
dist,ind = tree_ll.query([-300,300])
start_point = ll_points_cat[ind]
start_lonlat = ll_lonlat_cat[ind]
points_tmp.remove(start_point)
lonlat_tmp.remove(start_lonlat)
query_point = start_point
ll_points_merged.append(start_point)
ll_lonlat_merged.append(start_lonlat)

for i in range(len(ll_points_cat)-1):
    tree = KDTree(np.asarray(points_tmp))
    dist,ind = tree.query(start_point)
    
    dist2,ind2 = tree.query(query_point)
    
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




ll_points_merged = np.asarray(ll_points_merged, dtype=float)
ll_points_sorted = np.asarray(ll_points_sorted, dtype=float)
ll_points_smooth = gaussian_filter1d(ll_points_sorted, sigma=2/3, axis=0)
ll_points = ll_points_smooth

ll_lonlat_merged = np.asarray(ll_lonlat_merged, dtype=float)
ll_lonlat_sorted = np.asarray(ll_lonlat_sorted, dtype=float)
ll_lonlat_smooth = gaussian_filter1d(ll_lonlat_sorted, sigma=2/3, axis=0)
ll_lonlat = ll_lonlat_smooth




fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
plot_polygon(polygon, ax=ax, add_points=False, facecolor=None, edgecolor='k', linewidth=1)
ax.plot(x_cl, y_cl, 'k', linewidth=1.5)
ax.plot(ll_points_merged[:,0], ll_points_merged[:,1], 'r', linewidth=1)
ax.plot(ll_points_sorted[:,0], ll_points_sorted[:,1], 'b', linewidth=1)
ax.plot(ll_points_smooth[:,0], ll_points_smooth[:,1], 'w', linewidth=1)
# ax.scatter(zonal_points_filtered[:,0], zonal_points_filtered[:,1], marker='.', s=2, c='r')
# ax.scatter(merid_points_filtered[:,0], merid_points_filtered[:,1], marker='.', s=2, c='b')
# ax.plot(np.column_stack((cl_points_matched_zonal[:,0], x_zonal)).transpose(), np.column_stack((cl_points_matched_zonal[:,1], y_zonal)).transpose(), '-r', linewidth=0.75)
ax.set_title(f"QLCS polygon with centerline and zonal points")
ax.set_xlim([-200,100])
ax.set_ylim([-100,200])
plt.show()




if True:
    ll_zonal = {'x':x_ll_u, 'y':y_ll_u, 'i':idx_ll_u, 'j':jdy_ll_u,
                'x_filtered':zonal_points_filtered[:,0], 'y_filtered':zonal_points_filtered[:,1]}
    ll_meridional = {'x':x_ll_v, 'y':y_ll_v, 'i':idx_ll_v, 'j':jdy_ll_v,
                     'x_filtered':merid_points_filtered[:,0], 'y_filtered':merid_points_filtered[:,1]}
    data = {'ll_points':ll_points_smooth, 'll_lonlat':ll_lonlat_smooth, 'zonal':ll_zonal, 'meridional':ll_meridional, 'u_sm':u_sm, 'v_sm':v_sm}
    
    dbfile = open(fp+'leading_line_test_2.pkl', 'wb')
    pickle.dump(data, dbfile)
    dbfile.close()




#%% Get inflow polygon

dbfile = open(fp+'leading_line_test_2.pkl', 'rb')
tmp = pickle.load(dbfile)
ll_points = tmp['ll_points']
x_ll = ll_points[:,0]
y_ll = ll_points[:,1]
dbfile.close()

theta_local = np.zeros((len(ll_points),), dtype=float)
theta_norm = np.zeros((len(ll_points),), dtype=float)
x_proj = np.zeros((len(ll_points),), dtype=float)
y_proj = np.zeros((len(ll_points),), dtype=float)
dist = 40

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




for i in range(len(theta_local)):
    if abs(u_sm) > abs(v_sm):
        if u_sm > 0:
            theta_norm[i] = theta_local[i] + np.pi/2
        elif u_sm < 0:
            theta_norm[i] = theta_local[i] - np.pi/2
    elif abs(v_sm) > abs(u_sm):
        if v_sm > 0:
            theta_norm[i] = theta_local[i] + np.pi/2
        elif v_sm < 0:
            theta_norm[i] = theta_local[i] - np.pi/2
    
    
    x_proj[i] = x_ll[i] + dist*np.cos(theta_norm[i])
    y_proj[i] = y_ll[i] + dist*np.sin(theta_norm[i])



# QC polygon

qc_points = ll_points.tolist()
qc_points.append([-300,-300])
qc_polygon = Polygon(qc_points)


# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
# plot_polygon(qc_polygon, ax=ax, add_points=False, facecolor=None, edgecolor='k', linewidth=1)
# ax.plot(x_ll, y_ll, 'k', linewidth=1)
# ax.plot(x_proj, y_proj, 'k', linewidth=1)
# ax.set_title(f"Projection QC polygon")
# ax.set_xlim([-300,300])
# ax.set_ylim([-300,300])
# plt.show()


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



# ll_points = ll_points.tolist()
leading_points = [[x_ll[i], y_ll[i]] for i in range(len(ll_points))]
proj_points = [[x_proj_qc[i], y_proj_qc[i]] for i in range(len(ll_points))][::-1]

inflow_polygon = shapely.buffer(Polygon(leading_points+proj_points), 0)


fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
plot_polygon(inflow_polygon, ax=ax, add_points=False, facecolor=None, edgecolor='k', linewidth=1)
ax.set_title(f"Inflow polygon")
ax.set_xlim([-300,300])
ax.set_ylim([-300,300])
plt.show()



#%% Inflow shear analysis

# need to save lonm and latm to pickle probably. need those for era5
fp = 'C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/tor_outbreaks/'

# yyyyt = 2026
# mmt = 8
# ddt = 2
hht = int(timestr[:2])
# timestr = '161342'
timt = f"{yyyyt}-{mmt:02.0f}-{ddt:02.0f}T{hht:02.0f}:00:00.000000000"

fn_preslev = fp + f"era5_{yyyyt}{mmt}{ddt}_preslevs.nc"
fn_singlev = fp + f"era5_{yyyyt}{mmt}{ddt}_singlevs.nc"

datap = xr.open_dataset(fn_preslev)
datas = xr.open_dataset(fn_singlev)

latitude = datap['latitude'][:].values
longitude = datap['longitude'][:].values

lati = slice(np.argmin(abs(latitude-np.min(glat))), np.argmin(abs(latitude-np.max(glat)))+1)
loni = slice(np.argmin(abs(longitude-np.min(glon))), np.argmin(abs(longitude-np.max(glon)))+1)
latt = latitude[lati]
lont = longitude[loni]

data01 = datap.sel(latitude=slice(latt[0],latt[-1]), longitude=slice(lont[0],lont[-1]), valid_time=timt)
data02 = datas.sel(latitude=slice(latt[0],latt[-1]), longitude=slice(lont[0],lont[-1]), valid_time=timt)

datap.close()
datas.close()

z = data01['z'].values/9.81
u = data01['u'].values
v = data01['v'].values
orog = data02['z'].values/9.81
u10 = data02['u10'].values
v10 = data02['v10'].values

data01.close()
data02.close()



for j in range(len(latt)):
    for i in range(len(lont)):
        point = Point()
    point = Point(x_proj[i], y_proj[i])
    is_inside = qc_polygon.contains(point)


# ERA5 31-km resolution might be a problem here. might need to interpolate onto a finer grid idkkkkkkkk

































#%% empty space here

























#%% QLCS object ID testing and debugging


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
ax.set_title('Gridded CRef, prev')
plt.show()

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
ax.set_title('Gridded CRef')
plt.show()

fig,ax = plt.subplots(1, 1, figsize=(8,6))
plot_cfill(xm, ym, cref_smooth_next, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
ax.set_title('Gridded CRef, next')
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

#%%
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

#%%

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





#%% 





















