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
# latlon2xy
# get_storm_motion
# get_centroid
# get_leading_line
# get_inflow_polygon
# get_shear
# InterpolateToHeightAboveGround
# get_line_normal_and_parallel_shear
# get_local_tortuosity
# calc_QTor

# read_nexrad
# read_era5


### Functions to write
# read_eccc
# read_hrdps




#%% Load data from target and previous volumes + Cressman filtering of composite reflectivity

fp = 'C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/qlcs_tornado_outbreaks/'

radar_name = 'KBUF'
volume_time = datetime(2026, 8, 2, 16, 28, 24)

max_rng = 300 #set max range, km

# Load data files

if radar_name[0] == 'K':
    # Analysis volume
    filename = fp + f"{radar_name}/{radar_name}{yyyyt}{mmt:02.0f}{ddt:02.0f}_{timestr}_V06.ar2v"
    print(f"...Reading {filename}")
    cref, time, radar_lat, radar_lon, gx, gy, glat, glon = read_nexrad(filename, max_rng=max_rng)
    
    # Previous volume for storm motion estimation
    filenames = glob(fp+f"{radar_name}/{radar_name}{yyyyt}{mmt:02.0f}{ddt:02.0f}_*_V06.ar2v")
    ind = [i for i in range(len(filenames)) if timestr in filenames[i]][0]
    print(f"...Reading {filenames[ind-1]}")
    cref_prev, time_prev, _, _, _, _, _, _ = read_nexrad(filenames[ind-1], max_rng=max_rng)

elif radar_name[0] == 'C':
    df = pd.read_csv('C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/ECCC_radar_locations.csv',
                     sep=",", header=0, usecols=["Call sign", "Latitude", "Longitude"], index_col="Call sign")
    radar_lat = df.loc[radar_name]['Latitude']
    radar_lon = df.loc[radar_name]['Longitude']



#%%

# Cressman filtering
hres = 3 #gridded resolution, km
radius = 3 #Cressman radius of influence, km

# gx = gate_x[:,(rng<=max_rng)].astype(np.float32) #limit gates to range <= max range
# gy = gate_y[:,(rng<=max_rng)].astype(np.float32)
# glat = gate_lat[:,(rng<=max_rng)].astype(np.float32)
# glon = gate_lon[:,(rng<=max_rng)].astype(np.float32)
## xy and lat/lon meshgrids
xm,ym = np.meshgrid(np.arange(-max_rng, max_rng+hres, hres), np.arange(-max_rng, max_rng+hres, radius))
lonm,latm = np.meshgrid(np.linspace(np.min(glon), np.max(glon), len(xm)), np.linspace(np.min(glat), np.max(glat), len(xm))) #gate lat/lon mesh

print(f"...Cressman filtering analysis volume")
cref_smooth = cressman_interpolation_dask(gx, gy, cref, xm, ym, radius, chunk_size=5000)

print(f"...Cressman filtering previous volume")
cref_smooth_prev = cressman_interpolation_dask(gx, gy, cref_prev, xm, ym, radius, chunk_size=5000)



# # Originally used Gaussian interpolation until I figured out Dask
# pts = np.array([gate_x[:,(rng<=max_rng)].ravel(), gate_y[:,(rng<=max_rng)].ravel()]).transpose()
# vals = cref[:,(rng<=max_rng)].ravel().data
# xm,ym = np.meshgrid(np.arange(-max_rng, max_rng+hres, hres), np.arange(-max_rng, max_rng+hres, hres))
# xi = np.array([xm.ravel(), ym.ravel()]).transpose()
# cref_interp_flat = interpolate_to_points(pts, vals, xi, interp_type='linear', search_radius=hres)
# cref_interp = cref_interp_flat.reshape(xm.shape)
# cref_smooth = filters.gaussian(cref_interp, sigma=0.8)


# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
# ax.set_title('Cressman interpolated cref')
# plt.show()


#%% Retrieve QLCS object ID for analysis and previous volumes


qlcs_labels, qlcs_regions = get_qlcs_objects(cref_smooth, hres)
qlcs_labels_prev, qlcs_regions_prev = get_qlcs_objects(cref_smooth_prev, hres)

#% if more than one object in the analysis volume, pick the one with the highest cref - this may not be a permanent solution - maybe largest area?
if len(qlcs_regions) > 1:
    ## max cref
    crefmax = np.array([np.nanmax(cref_smooth[(qlcs_labels==qlcs_regions[i].label)]) for i in range(len(qlcs_regions))])
    qlcs_region = qlcs_regions[np.argmax(crefmax)]
    ## max area
    # areamax = np.array([qlcs_regions[i].area for i in range(len(qlcs_regions))])
    # qlcs_region = qlcs_regions[np.argmax(areamax)]
    
    qlcs_obj = np.where(qlcs_labels == qlcs_region.label, 1, 0)
else:
    qlcs_obj = np.where(qlcs_labels>0, 1, 0)
    qlcs_region = qlcs_regions[0]



# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, qlcs_labels, 'dbz', ax, datalims=[0,len(qlcs_regions)], cmap='HomeyerRainbow')
# ax.set_title('Final QLCS objects\n filtered for max dBZ, merged, and filtered for length and eccentricity ')
# plt.show()

# fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(xm, ym, qlcs_obj, 'dbz', ax, datalims=[0,1], cmap='HomeyerRainbow')
# ax.set_title('Final QLCS object\n Object with highest dBZ ')
# plt.show()




#%% Get centroids from previous and next scan and estimate storm motion

if True:
    # Get centroids
    centroid = get_centroid(qlcs_region, xm, ym, reference_centroid=None)
    centroid_prev = get_centroid(qlcs_regions_prev, xm, ym, reference_centroid=centroid)
    x_centroid, y_centroid = centroid[0], centroid[1]
    x_centroid_prev, y_centroid_prev = centroid_prev[0], centroid_prev[1]
    
    [u_sm, v_sm] = get_storm_motion(centroid, centroid_prev, time, time_prev)
    
if False:
    # Analysis volume
    j_centroid,i_centroid = qlcs_region.centroid
    jc1,ic1,jc2,ic2 = np.floor(j_centroid), np.floor(i_centroid), np.ceil(j_centroid), np.ceil(i_centroid)
    j1,i1,j2,i2 = int(jc1), int(ic1), int(jc2), int(ic2)
    x_centroid = (1 - abs(i_centroid-ic1))*xm[j1,i1] + (1 - abs(i_centroid-ic2))*xm[j1,i2]
    y_centroid = (1 - abs(j_centroid-jc1))*ym[j1,i1] + (1 - abs(j_centroid-jc2))*ym[j2,i1]
    centroid = [x_centroid, y_centroid]
    
    # Previous volume
    qlcs_obj_prev = np.zeros(shape=qlcs_labels_prev.shape, dtype=int)
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
        centroid_prev = [x_centroid_prev, y_centroid_prev]
    else:
        if isinstance(region, list):
            j_centroid,i_centroid = qlcs_regions_prev[0].centroid
        else:
            j_centroid,i_centroid = qlcs_regions_prev.centroid
        jc1,ic1,jc2,ic2 = np.floor(j_centroid), np.floor(i_centroid), np.ceil(j_centroid), np.ceil(i_centroid)
        j1,i1,j2,i2 = int(jc1), int(ic1), int(jc2), int(ic2)
        x_centroid_prev = (1 - abs(i_centroid-ic1))*xm[j1,i1] + (1 - abs(i_centroid-ic2))*xm[j1,i2]
        y_centroid_prev = (1 - abs(j_centroid-jc1))*ym[j1,i1] + (1 - abs(j_centroid-jc2))*ym[j2,i1]
        centroid_prev = [x_centroid_prev, y_centroid_prev]
        qlcs_obj_prev = qlcs_labels_prev
    
    # Storm motion
    dt_prev = (time - time_prev).total_seconds()
    u_sm = (x_centroid - x_centroid_prev)*1000 / dt_prev
    v_sm = (y_centroid - y_centroid_prev)*1000 / dt_prev




#%% Get leading line and inflow polygon

ll_points, ll_lonlat = get_leading_line(qlcs_obj, [u_sm, v_sm], xm, ym, latm, lonm)

inflow_polygon, theta_norm, theta_local = get_inflow_polygon(ll_points, [u_sm, v_sm])


#%% Load ERA5 data for inflow shear analysis

fp2 = 'C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/tor_outbreaks/'


latlims = [np.min(latm), np.max(latm)]
lonlims = [np.min(lonm), np.max(lonm)]
data, latt, lont = read_era5(time, fp2, latlims, lonlims)

p = data['p']
z = data['z']
orog = data['orog']
u = data['u']
v = data['v']
u10 = data['u10']
v10 = data['v10']


lonm2,latm2 = np.meshgrid(lont, latt)
x = np.zeros(lonm2.shape, dtype=float)
y = np.zeros(lonm2.shape, dtype=float)

for j in range(len(long)):
    xtmp,ytmp = latlon2xy(latm2[j,:], lonm2[j,:], radar_lat, radar_lon)
    x[j,:] = xtmp
    y[j,:] = ytmp



#%%
xm2,ym2 = np.meshgrid(np.arange(-max_rng, max_rng+hres, hres), np.arange(-max_rng, max_rng+hres, hres))

u_interp = np.zeros(shape=(len(p),xm.shape[1],xm.shape[0]), dtype=float)
v_interp = np.zeros(shape=(len(p),xm.shape[1],xm.shape[0]), dtype=float)
z_interp = np.zeros(shape=(len(p),xm.shape[1],xm.shape[0]), dtype=float)

for k in range(len(p)):
    print(f"Pressure level {p[k]:.0f} hPa")
    u_interp[k,:,:] = cressman_interpolation_dask(x, y, u[k,:,:], xm, ym, 62, chunk_size=3000)
    v_interp[k,:,:] = cressman_interpolation_dask(x, y, v[k,:,:], xm, ym, 62, chunk_size=3000)
    z_interp[k,:,:] = cressman_interpolation_dask(x, y, z[k,:,:], xm, ym, 62, chunk_size=3000)
u10_interp = cressman_interpolation_dask(x, y, u10, xm, ym, 62, chunk_size=3000)
v10_interp = cressman_interpolation_dask(x, y, v10, xm, ym, 62, chunk_size=3000)
orog_interp = cressman_interpolation_dask(x, y, orog, xm, ym, 62, chunk_size=3000)


#%% Calculate QTor!!!

shear03 = get_shear(z_interp, orog_interp, u_interp, v_interp, u10_interp, v10_interp, 3000)
shear01 = get_shear(z_interp, orog_interp, u_interp, v_interp, u10_interp, v10_interp, 1000)



points = np.column_stack( (xm.ravel(), ym.ravel()))


LN03, LP01 = get_line_normal_and_parallel_shear(shear03, shear01, ll_points, theta_norm, inflow_polygon, points)
tortuosity = get_local_tortuosity(ll_points)
qtor = calc_QTor(LN03, LP01, tortuosity)



#%%

from scipy.interpolate import griddata

qtor_grid = griddata(ll_points, qtor, (xm,ym), fill_value=0, method='cubic')

fig,ax = plt.subplots(1, 1, figsize=(8,6))
# plot_cfill(gx, gy, cref, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
plot_cfill(xm, ym, cref_smooth, 'dbz', ax, datalims=[0,70], cmap='HomeyerRainbow')
ax.contour(xm, ym, qtor_grid, levels=[0.5,1.0], colors=['k'], linewidths=[0.5,1])
ax.set_title('Original composite reflectivity + QTor')
plt.show()




#%% empty space here



















#%% Shear/tortuosity analysis and calculate QTor - development

shear03_vals = np.column_stack( (shear03[0].ravel(), shear03[1].ravel()))
shear01_vals = np.column_stack( (shear01[0].ravel(), shear01[1].ravel()))

tree = KDTree(ll_points)

mask = [False] * len(points)
for i in range(len(points)):
    point = Point(points[i])
    mask[i] = inflow_polygon.covers(point)

inflow_points = points[mask]
inflow_shear03 = shear03_vals[mask]
inflow_shear01 = shear01_vals[mask]

dists,inds = tree.query(inflow_points) # Find nearest leading line point in tree to each inflow point


LN03 = np.zeros((len(ll_points),), dtype=float)
LP01 = np.zeros((len(ll_points),), dtype=float)
tortuosity = np.zeros((len(ll_points),), dtype=float)
# t_factor = np.zeros((len(ll_points),), dtype=float)
# qtor = np.zeros((len(ll_points),), dtype=float)

for i in range(len(ll_points)):
    if np.any(inds == i):
        s03 = inflow_shear03[(inds==i)]
        s01 = inflow_shear01[(inds==i)]
        S03 = np.sqrt(s03[:,0]**2 + s03[:,1]**2)
        S01 = np.sqrt(s01[:,0]**2 + s01[:,1]**2)
        
        orientation_local = [np.cos(theta_norm[i]+np.pi/2), np.sin(theta_norm[i]+np.pi/2)] #line orientation to the left of normal vector
        
        lp01 = np.zeros((len(s03),), dtype=float)
        ln03 = np.zeros((len(s03),), dtype=float)
        for j in range(len(s03)):
            s03_norm = s03[j] / S03[j]
            s01_norm = s01[j] / S01[j]
            lp03_percent = np.abs(np.dot(orientation_local, s03_norm))
            lp01_percent = np.abs(np.dot(orientation_local, s01_norm))
            
            lp01[j] = lp01_percent * S01[j]
            ln03_mag = (1 - lp03_percent) * S03[j]
            ln03_sgn = np.sign(np.cross(s03_norm, orientation_local))
            ln03[j] = ln03_mag * ln03_sgn
        
        LP01[i] = np.nanmean(lp01)
        LN03[i] = np.nanmean(ln03)
    
        
    # Local tortuosity - moving block
    
    inds_left = np.arange(i-5,i+1)
    inds_right = np.arange(i,i+6)
    
    ileft = inds_left[(inds_left>=0)][0]
    iright = inds_right[(inds_right<len(ll_points))][-1]
    
    segment_length = 0.0
    for k in np.arange(ileft,iright):
        segment_length += distance(ll_points[k], ll_points[k+1])
    
    endpoint_length = distance(ll_points[ileft], ll_points[iright])
    
    tortuosity[i] = segment_length / endpoint_length

LN_factor = LN03 / 10.0
LP_factor = LP01 / 8.0
LN_factor[(LN_factor<0.5)] = 0.5 #minimum of 0.5 for shear terms
LP_factor[(LP_factor<0.5)] = 0.5
LN_factor[(LN_factor>3.0)] = 3.0 #shear terms capped at 3.0
LP_factor[(LP_factor>3.0)] = 3.0
T_factor = np.where(tortuosity > 1.05, 2, 1)


# Calculate QTor!

qtor = LN_factor * LP_factor * T_factor







#%% Initial leading line identification using storm motion -- development



# E-W search
x_zonal = []
y_zonal = []
lat_zonal = []
lon_zonal = []

for jdy in range(xm.shape[0]):
    if np.any(qlcs_obj[jdy,:] > 0):
        if u_sm > 0:
            idx = np.where(qlcs_obj[jdy,:]>0)[0][-1]
        elif u_sm < 0:
            idx = np.where(qlcs_obj[jdy,:]>0)[0][0]
        
        x_zonal.append(xm[jdy,idx])
        y_zonal.append(ym[jdy,idx])
        lat_zonal.append(latm[jdy,idx])
        lon_zonal.append(lonm[jdy,idx])
        

# N-S search
x_merid = []
y_merid = []
lat_merid = []
lon_merid = []

for idx in range(ym.shape[1]):
    if np.any(qlcs_obj[:,idx] > 0):
        if v_sm > 0:
            jdy = np.where(qlcs_obj[:,idx]>0)[0][-1]
        elif v_sm < 0:
            jdy = np.where(qlcs_obj[:,idx]>0)[0][0]
        
        x_merid.append(xm[jdy,idx])
        y_merid.append(ym[jdy,idx])
        lat_merid.append(latm[jdy,idx])
        lon_merid.append(lonm[jdy,idx])


#% Find centerline
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



if False:
    data_cl = {'x':x_cl, 'y':y_cl, 'centerline':centerline, 'trimmed_centerline':trimmed_centerline}
    dbfile = open(fp+'centerline_test.pkl', 'wb')
    pickle.dump(data_cl, dbfile)
    dbfile.close()



#% QC and smooth leading line points
#   Find closest centerline point to each leading line point,
#   create vector between centerline and leading line points,
#   exclude any leading line points misaligned with expected SM
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







#%% Get inflow polygon -- development

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





















