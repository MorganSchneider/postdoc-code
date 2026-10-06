# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 12:41:29 2026

@author: mschne28
"""

from QTORutils import *


#%% 


fp = 'C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/qlcs_tornado_outbreaks/'

radar_name = 'KDTX'
year = 2026
month = 9
day = 3
hour = 20

model_data = "ERA5" # ERA5, HRDPS, HRRR
model_run = 18 #UTC, for HRDPS, HRRR

use_hdf5 = True


datestr = f"{year}{month:02.0f}{day:02.0f}"

max_rng = 300 #set max range, km
hres = 3 #gridded resolution, km
radius = 3 #Cressman radius of influence, km
xm,ym = np.meshgrid(np.arange(-max_rng, max_rng+hres, hres), np.arange(-max_rng, max_rng+hres, radius))


# ECCC
if radar_name[0] == 'C':
    if use_hdf5:
        files = glob(fp + f"{radar_name}/{datestr}T{hour:02.0f}*_{radar_name}.hdf5")
    else:
        files = glob(fp + f"{radar_name}/*.geojson")
# NEXRAD
elif radar_name[0] == 'K':
    files = glob(fp + f"{radar_name}/{radar_name}{datestr}_{hour:02.0f}*_V06.ar2v")


for file in files:
    
    # ECCC
    if radar_name == 'C':
        if use_hdf5:
            # Analysis volume
            volume_time = datetime.strptime(os.path.basename(file), f"%Y%m%dT%H%MZ_MSC_Radar-VolumeScans_{radar_name}.hdf5")
            timestr = volume_time.strftime("%H%M")
            
            print(f"...Reading {os.path.basename(file)}")
            cref_smooth, time, radar_lat, radar_lon, gx, gy, glat, glon = read_eccc_hdf5(file, xm, ym, max_rng=max_rng, radius=radius)
            
            # Previous volume for storm motion estimation
            files_all = glob(fp+f"{radar_name}/{datestr}T*_{radar_name}.hdf5")
            ind = [i for i in range(len(files_all)) if timestr in files_all[i]][0]
            print(f"...Reading {os.path.basename(files_all[ind-1])}")
            cref_smooth_prev, time_prev, _, _, _, _, _, _ = read_eccc_hdf5(files_all[ind-1], xm, ym, max_rng=max_rng, radius=radius)
        else:
            volume_time = datetime.strptime(os.path.basename(file), f"%Y%m%dT%H%MZ_MSC_Radar-VolumeScans_{radar_name}.geojson")
    
    # NEXRAD
    elif radar_name[0] == 'K':
        # Analysis volume
        volume_time = datetime.strptime(os.path.basename(file), f"{radar_name}%Y%m%d_%H%M%S_V06.ar2v")
        timestr = volume_time.strftime("%H%M%S")
        
        print(f"...Reading {os.path.basename(file)}")
        cref, time, radar_lat, radar_lon, gx, gy, glat, glon = read_nexrad(file, max_rng=max_rng)
        
        # Previous volume for storm motion estimation
        files_all = glob(fp+f"{radar_name}/{radar_name}{datestr}_*_V06.ar2v")
        ind = [i for i in range(len(files_all)) if timestr in files_all[i]][0]
        print(f"...Reading {os.path.basename(files_all[ind-1])}")
        cref_prev, time_prev, _, _, _, _, _, _ = read_nexrad(files_all[ind-1], max_rng=max_rng)
        
        print(f"...Cressman filtering analysis volume")
        cref_smooth = cressman_interpolation_dask(gx, gy, cref, xm, ym, radius, chunk_size=5000)

        print(f"...Cressman filtering previous volume")
        cref_smooth_prev = cressman_interpolation_dask(gx, gy, cref_prev, xm, ym, radius, chunk_size=5000)
    
    
    lonm,latm = np.meshgrid(np.linspace(np.min(glon), np.max(glon), len(xm)), np.linspace(np.min(glat), np.max(glat), len(xm))) #gate lat/lon mesh

    
    print(f"...Identifying QLCS objects")
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
    
    
    print("...Estimating storm motion")
    centroid = get_centroid(qlcs_region, xm, ym, reference_centroid=None)
    centroid_prev = get_centroid(qlcs_regions_prev, xm, ym, reference_centroid=centroid)
    x_centroid, y_centroid = centroid[0], centroid[1]
    x_centroid_prev, y_centroid_prev = centroid_prev[0], centroid_prev[1]
    [u_sm, v_sm] = get_storm_motion(centroid, centroid_prev, time, time_prev)
    
    print("...Finding leading line and inflow polygon")
    ll_points, ll_lonlat, obj_contour = get_leading_line(qlcs_obj, [u_sm, v_sm], xm, ym, latm, lonm)
    inflow_polygon, theta_norm, theta_local = get_inflow_polygon(ll_points, [u_sm, v_sm])
    
    # Get wind field data
    if file == files[0]:
        if model_data == "ERA5":
            print("...Reading ERA5 data")
            fp2 = 'C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/tor_outbreaks/'
            latlims = [np.min(latm), np.max(latm)]
            lonlims = [np.min(lonm), np.max(lonm)]
            rad = 62
            
            data, latt, lont = read_era5_netcdf(time, fp2, latlims, lonlims)
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
            for j in range(len(lonm2)):
                xtmp,ytmp = latlon2xy(latm2[j,:], lonm2[j,:], radar_lat, radar_lon)
                x[j,:] = xtmp
                y[j,:] = ytmp
            
            # u = np.zeros(shape=(len(p),xm.shape[1],xm.shape[0]), dtype=float)
            # v = np.zeros(shape=(len(p),xm.shape[1],xm.shape[0]), dtype=float)
            # z = np.zeros(shape=(len(p),xm.shape[1],xm.shape[0]), dtype=float)
            # for k in range(len(p_era5)):
            #     print(f"Pressure level {p[k]:.0f} hPa")
            #     u[k,:,:] = cressman_interpolation_dask(x, y, u[k,:,:], xm, ym, 62, chunk_size=3000)
            #     v[k,:,:] = cressman_interpolation_dask(x, y, v[k,:,:], xm, ym, 62, chunk_size=3000)
            #     z[k,:,:] = cressman_interpolation_dask(x, y, z[k,:,:], xm, ym, 62, chunk_size=3000)
            # u10 = cressman_interpolation_dask(x, y, u10, xm, ym, 62, chunk_size=3000)
            # v10 = cressman_interpolation_dask(x, y, v10, xm, ym, 62, chunk_size=3000)
            # orog = cressman_interpolation_dask(x, y, orog, xm, ym, 62, chunk_size=3000)
        
        elif model_data == "HRDPS":
            fcst_hour = time.hour - model_run
            print(f"...Reading {model_run:02.0f}Z T{fcst_hour:03.0f} HRDPS data")
            latlims = [np.min(latm), np.max(latm)]
            lonlims = [np.min(lonm), np.max(lonm)]
            rad = radius
    
            fp2 = f"C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/qlcs_tornado_outbreaks/hrdps/{datestr}/{model_run:02.0f}z/T{fcst_hour:03.0f}/"
    
            data, latt, lont = read_hrdps(time, model_run, fcst_hour, fp2, latlims, lonlims)
            p = data['p']
            z = data['z']
            orog = data['orog']
            u = data['u']
            v = data['v']
            u10 = data['u10']
            v10 = data['v10']
    
            x = np.zeros(lont.shape, dtype=float)
            y = np.zeros(lont.shape, dtype=float)
    
            for j in range(lont.shape[0]):
                xtmp,ytmp = latlon2xy(latt[j,:], lont[j,:], radar_lat, radar_lon)
                x[j,:] = xtmp
                y[j,:] = ytmp

        elif model_data == "HRRR":
            fcst_hour = time.hour - model_run
            print(f"...Reading {model_run:02.0f}Z T{fcst_hour:03.0f} HRRR data")
            latlims = [np.min(latm), np.max(latm)]
            lonlims = [np.min(lonm), np.max(lonm)]
            rad = radius

            fp2 = f"C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/qlcs_tornado_outbreaks/hrrr/{datestr}/{model_run:02.0f}z/"
            
            data, latt, lont = read_hrrr(time, model_run, fcst_hour, fp2, latlims, lonlims)
            p = data['p']
            z = data['z']
            orog = data['orog']
            u = data['u']
            v = data['v']
            u10 = data['u10']
            v10 = data['v10']
                
            x = np.zeros(lont.shape, dtype=float)
            y = np.zeros(lont.shape, dtype=float)
                
            for j in range(lont.shape[0]):
                xtmp,ytmp = latlon2xy(latt[j,:], lont[j,:], radar_lat, radar_lon)
                x[j,:] = xtmp
                y[j,:] = ytmp
    
    
    print("...Calculating wind shear")
    shear03 = get_shear(z, orog, u, v, u10, v10, 3000)
    shear01 = get_shear(z, orog, u, v, u10, v10, 1000)
    ushear03 = cressman_interpolation_dask(x, y, shear03[0], xm, ym, rad, chunk_size=3000)
    vshear03 = cressman_interpolation_dask(x, y, shear03[1], xm, ym, rad, chunk_size=3000)
    ushear01 = cressman_interpolation_dask(x, y, shear01[0], xm, ym, rad, chunk_size=3000)
    vshear01 = cressman_interpolation_dask(x, y, shear01[1], xm, ym, rad, chunk_size=3000)
    shear03 = (ushear03, vshear03)
    shear01 = (ushear01, vshear01)
    
    print("...Calculating QTor")
    points = np.column_stack( (xm.ravel(), ym.ravel()))
    LN03, LP01 = get_line_normal_and_parallel_shear(shear03, shear01, ll_points, theta_norm, inflow_polygon, points)
    tortuosity = get_local_tortuosity(ll_points)
    qtor = calc_QTor(LN03, LP01, tortuosity)
    qtor_grid = griddata(ll_points, qtor, (xm,ym), fill_value=0, method='cubic')
    
    save_data = dict(qtor=qtor, ln03=LN03, lp01=LP01, tortuosity=tortuosity, xm=xm, ym=ym, ll_points=ll_points, cref_smooth=cref_smooth, qlcs_obj=qlcs_obj)
    
    dbfile = open(fp + f"qtor_{datestr}_{timestr}.pkl", 'wb')
    pickle.dump(dbfile, save_data)
    dbfile.close()

































