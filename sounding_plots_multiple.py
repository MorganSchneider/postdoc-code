# -*- coding: utf-8 -*-
"""
Created on Fri May  1 11:47:09 2026

@author: L. Schielicke
"""

### Import packages ###

from era5utils import *


#%% Choose time and location range for sounding

#Didsbury, AB
#Date: 7/1/2023
#Coordinates: (51.6107,-114.198)
#Time (UTC): 17:45

yyyyt = 2023
mmt   = 7
ddt   = 1
hht   = 20
lattstart=51.6107  # from NTP/NHP website
lontstart=-114.198 # from NTP/NHP website
case = 'tornado' #'hail'



inversionfile = 0   # set inversionfile = 0 if you do not need sounding files with a low-level inversion, 
                    # set inversionfile to any other value if you want to write out inversion profiles - check values!
time_str= "%d%s%sT%s" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))
time_str2= "%d%s%s_%s" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))
pngfolder = "./input_soundings/%d%s%s_%s/" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))
#geojsonfolder = pngfolder+"geojson/" 
geojsonfolder = pngfolder+"geojson_" 
os.makedirs(pngfolder, exist_ok=True)
#os.makedirs(geojsonfolder, exist_ok=True)

data_file_preslev=f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_preslevs.nc"
data_file_singlev=f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_singlevs.nc"

name='Didsbury'#'Alonsa_MB'



# githubfolder=f"https://raw.githubusercontent.com/severe-weather-lab/soundings/refs/heads/main/{time_str}/" # OLD!!!
basefolder="https://raw.githubusercontent.com/severe-weather-lab/soundings/refs/heads/main/"
githubfolder=f"{basefolder}input_multiple_soundings/{case}/{name}/{time_str2}/"
print(githubfolder)



latt_south = np.round(lattstart/0.25)*0.25-2.5 #49.5       # southern latitude of the domain
lont_west  = np.round(lontstart/0.25)*0.25-2.5     # western longitude of the domain 
print(latt_south,lont_west)

timt="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))
#timt="2023-07-01T17:00:00.000000000"
print(timt)



#%% Open dataset

data = xr.open_dataset(data_file_preslev)
datas = xr.open_dataset(data_file_singlev)

data
datas

#%% Determine closest lat/lon of event for sounding plots
# Fields plotted +/- 2.5 deg around location of event

latitude = data['latitude'][:].values
longitude= data['longitude'][:].values

loni=np.argmin(np.abs(longitude-lontstart))
lati=np.argmin(np.abs(latitude-lattstart))
print(lati,loni)
print(latitude[lati],longitude[loni])



#%% Loop over lat/lon range and save sounding plots

for i in range(0,21):
    lont=longitude[loni]-2.5+i*0.25
    print(lont)
    for j in range(0,21):
        latt=latitude[lati]-2.5+j*0.25
        p,z,T,q,theta,Td,u,v,u10,v10,speed,direc,cape,cin,sfcp,orog,q2m,theta2m,td2m,t2m,leftm,meanm,rightm,parcel_prof,lcl_pressure,lcl_temperature = extract_data(latt,lont,timt,data,datas)
        
        pngfolder=None
        plot_with_hodograph(pngfolder,timt,lont,latt,z,orog,T,Td,p,u,v,u10,v10,leftm,meanm,rightm,lcl_pressure, lcl_temperature,parcel_prof,pngfolder)


#%% Interpolate and plot

u500,v500,q500,lapse_rate= PlotMoistureMap1(data,datas,longitude[loni]-2.5,longitude[loni]+2.5,latitude[lati]-2.5,latitude[lati]+2.5,case,1)



#%% Write GeoJSON files for ArcGIS interactive map

# Wind 500 m AGL -- vectors, saved as LineString
features = []

for i in range(0,21):
    lont=longitude[loni]-2.5+i*0.25
    for j in range(0,21):
        latt=latitude[lati]-2.5+j*0.25
        u_val = u500[lati+10-j, loni-10+i]
        v_val = v500[lati+10-j, loni-10+i]
        
        # Skip NaNs or calm wind
        if np.isnan(u_val) or np.isnan(v_val):
            continue

        # Convert meters to degrees
        deg_per_meter_lat = 1 / 111000
        deg_per_meter_lon = 1 / (111000 * np.cos(np.radians(latt)))

        # Arrow scale factor (adjust as needed)
        arrow_scale = 5000

        dx = u_val * arrow_scale * deg_per_meter_lon
        dy = v_val * arrow_scale * deg_per_meter_lat

        start = (lont, latt)
        end = (lont + dx, latt + dy)

        geometry = LineString([start, end])

        speed = np.sqrt(u_val**2 + v_val**2)
        direction = (180 / np.pi) * np.arctan2(u_val, v_val)  # meteorological convention: from

        features.append({
            'geometry': geometry,
            'u': u_val,
            'v': v_val,
            'speed': speed,
            'direction': direction,
            'timestamp': time_str
        })

gdf = gpd.GeoDataFrame(features, crs='EPSG:4326')
    
out_path = geojsonfolder+f'wind500m_{time_str}.geojson'
gdf.to_file(out_path, driver='GeoJSON')
print(f'Saved: {out_path}')




# Specific humidity 500 m AGL -- saved as Point
features = []

for i in range(0,21):
    lont=longitude[loni]-2.5+i*0.25
    for j in range(0,21):
        latt=latitude[lati]-2.5+j*0.25
        val = q500[lati+10-j, loni-10+i]*1000
        lapse = lapse_rate[lati+10-j, loni-10+i]

        if np.isnan(val):
            continue

        geometry = Point(lont, latt)

        features.append({
            'geometry': geometry,
            'humidity@500m': val,
            'lapserate': lapse,
            'timestamp': time_str
        })

gdf = gpd.GeoDataFrame(features, crs='EPSG:4326')
out_path = geojsonfolder + f'q500m_points_{time_str}.geojson'
gdf.to_file(out_path, driver='GeoJSON')
print(f'Saved: {out_path}')




# Lapse rate contours -- contours, saved as LineString
features = []

X,Y = np.meshgrid(longitude,latitude)
Xs=X[lati-10:lati+11,loni-10:loni+11]
Ys=Y[lati-10:lati+11,loni-10:loni+11]
field=lapse_rate[lati-10:lati+11,loni-10:loni+11]*1000

# Define levels and styling
levels = {
    8.0: {"color": "purple", "style": "solid", "width": 1},
    7.0: {"color": "red", "style": "solid",  "width": 1},
    6.5: {"color": "red", "style": "dashed",   "width": 1},
    6.0: {"color": "red", "style": "dotted", "width": 0.5},
    5.0: {"color": "gray", "style": "dotted", "width": 0.1},
}

fig,ax= plt.subplots()
    
for level, style in levels.items():
    cs = ax.contour(Xs, Ys, field, levels=[level], colors=style["color"], linewidths=style["width"])

#    for path in cs.collections[0].get_paths():
    for path in cs.get_paths():
        coords = path.vertices
        if len(coords) < 3:
            continue
        line = LineString(coords)
        features.append({
            "geometry": line,
            "level": level,
            "color": style["color"],
            "line_style": style["style"],
            "line_width": style["width"]
        })

# Write to GeoJSON
gdf = gpd.GeoDataFrame(features, crs="EPSG:4326")
out_path = geojsonfolder + f"lapse_rate_{time_str}.geojson" 
gdf.to_file(out_path, driver="GeoJSON")
print(f'Saved: {out_path}')




# Links to sounding -- saved as Point
features = []

for i in range(0,21):
    lont=longitude[loni]-2.5+i*0.25
    for j in range(0,21):
        latt=latitude[lati]-2.5+j*0.25
        img_url = githubfolder+"input_sounding_%s_%s_%s.png"  %(timt[:13],lont,latt)
        features.append({
            'geometry': Point(lont, latt),
            'lat': latt,
            'lon': lont,
            'image_url': img_url,
            'timestamp': time_str
        })

gdf = gpd.GeoDataFrame(features, crs='EPSG:4326')
out_path = geojsonfolder + f'sounding_points_{time_str}.geojson'
gdf.to_file(out_path, driver='GeoJSON')
print(f'Saved: {out_path}')




# CAPE -- field, saved as point
features = []

data02 = datas.sel(valid_time=timt)
cape = data02["cape"]

for i in range(0,21):
    lont=longitude[loni]-2.5+i*0.25
    for j in range(0,21):
        latt=latitude[lati]-2.5+j*0.25
        val = np.array(cape[lati+10-j, loni-10+i])*1.
        if np.isnan(val):
            continue

        geometry = Point(lont, latt)

        features.append({
            'geometry': geometry,
            'cape': val,
            'timestamp': time_str
        })

gdf = gpd.GeoDataFrame(features, crs='EPSG:4326')
out_path = geojsonfolder + f'cape_points_{time_str}.geojson'
gdf.to_file(out_path, driver='GeoJSON')
print(f'Saved: {out_path}')




# Tornado or hail location -- saved as Point
features = []

geometry = Point(lontstart, lattstart)

features.append({
    'geometry': geometry,
    'case': case,
    'timestamp': time_str
})

gdf = gpd.GeoDataFrame(features, crs='EPSG:4326')
out_path = geojsonfolder + f'{case}_point_{time_str}.geojson'
gdf.to_file(out_path, driver='GeoJSON')
print(f'Saved: {out_path}')




# Check if link is working
print(img_url)












