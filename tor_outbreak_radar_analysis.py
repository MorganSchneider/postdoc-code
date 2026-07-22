# -*- coding: utf-8 -*-
"""
Created on Fri May 29 11:42:14 2026

@author: mschne28
"""

from era5utils import *

#%%


dbfile = open("C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/tor_outbreaks/tornado_locs.pkl", 'rb')
locs = pickle.load(dbfile)
dbfile.close()

### August 11, 2021 - Tornadoes 1, 3, 5 have possible ZDR columns
tilt = [0.0, 0.3135, 0.4670, 1.0740, 1.6890]

latt = locs['20210811']['loc1']['lat']
lont = locs['20210811']['loc1']['lon']
zdr_column_height_21Aug11_tor1 = getBeamHeight(latt, lont, tilt[4], radar_id='CASMR')

latt = locs['20210811']['loc3']['lat']
lont = locs['20210811']['loc3']['lon']
zdr_column_height_21Aug11_tor3 = getBeamHeight(latt, lont, tilt[3], radar_id='CASMR')

latt = locs['20210811']['loc5']['lat']
lont = locs['20210811']['loc5']['lon']
zdr_column_height_21Aug11_tor5 = getBeamHeight(latt, lont, tilt[4], radar_id='CASMR')
# other locs: 46.2837, -83.1276 at 1654 EDT, tilt 3;   46.3819, -82.5893 at 1724 EDT, tilt 4

zdr_col_hgts_21Aug11 = np.array([zdr_column_height_21Aug11_tor1, zdr_column_height_21Aug11_tor3, zdr_column_height_21Aug11_tor5])

print("---ZDR column heights, 2021-08-11---")
print(f"Tornado 1: {zdr_column_height_21Aug11_tor1:.3f} km")
print(f"Tornado 3: {zdr_column_height_21Aug11_tor3:.3f} km")
print(f"Tornado 5: {zdr_column_height_21Aug11_tor5:.3f} km")
# print(f"Mean:      {np.mean(zdr_col_hgts_21Aug11):.3f} km")
print(" ")




### July 23, 2025 - All 6 tornadoes have possible ZDR columns
tilt = [0.0, 0.1100, 0.5190, 1.0110]

latt = locs['20250623']['loc1']['lat']
lont = locs['20250623']['loc1']['lon']
zdr_column_height_25Jun23_tor12 = getBeamHeight(latt, lont, tilt[2], radar_id='CASMA')

latt = locs['20250623']['loc3']['lat']
lont = locs['20250623']['loc3']['lon']
zdr_column_height_25Jun23_tor3 = getBeamHeight(latt, lont, tilt[3], radar_id='CASMA')
# other locs: 48.2342, -74.2484 at 1754 EDT, tilt 3
latt = 48.2342
lont = -74.2484
zdr_column_height_25Jun23_tor3_2 = getBeamHeight(latt, lont, tilt[3], radar_id='CASMA')

latt = locs['20250623']['loc4']['lat']
lont = locs['20250623']['loc4']['lon']
zdr_column_height_25Jun23_tor4 = getBeamHeight(latt, lont, tilt[3], radar_id='CASMA')

tilt = [0.0, 0.5135, 0.9090, 1.3180]

latt = locs['20250623']['loc5']['lat']
lont = locs['20250623']['loc5']['lon']
zdr_column_height_25Jun23_tor5 = getBeamHeight(latt, lont, tilt[3], radar_id='CASSF')

latt = locs['20250623']['loc6']['lat']
lont = locs['20250623']['loc6']['lon']
zdr_column_height_25Jun23_tor6 = getBeamHeight(latt, lont, tilt[3], radar_id='CASSF')

zdr_col_hgts_25Jun23 = np.array([zdr_column_height_25Jun23_tor12, zdr_column_height_25Jun23_tor3, zdr_column_height_25Jun23_tor3_2,
                                 zdr_column_height_25Jun23_tor4, zdr_column_height_25Jun23_tor5, zdr_column_height_25Jun23_tor6])

print("---ZDR column heights, 2025-06-23---")
print(f"Tornado 1/2: {zdr_column_height_25Jun23_tor12:.3f} km")
print(f"Tornado 3: {zdr_column_height_25Jun23_tor3:.3f} km (1830 EDT)")
print(f"           {zdr_column_height_25Jun23_tor3_2:.3f} km (1754 EDT)")
print(f"Tornado 4: {zdr_column_height_25Jun23_tor4:.3f} km")
print(f"Tornado 5: {zdr_column_height_25Jun23_tor5:.3f} km")
print(f"Tornado 6: {zdr_column_height_25Jun23_tor6:.3f} km")
# print(f"Mean:      {np.mean(zdr_col_hgts_25Jun23):.3f} km")
print(" ")




### May 30, 2022 - Tornadoes 1, 2, 3 have possible weak ZDR columns
tilt = [0.0, 0.4805, 0.8155, 1.2935, 1.6865]

latt = locs['20220530']['loc1']['lat']
lont = locs['20220530']['loc1']['lon']
zdr_column_height_22May30_tor1 = getBeamHeight(latt, lont, tilt[4], radar_id='CASDR')

latt = locs['20220530']['loc2']['lat']
lont = locs['20220530']['loc2']['lon']
zdr_column_height_22May30_tor2 = getBeamHeight(latt, lont, tilt[4], radar_id='CASDR')

latt = locs['20220530']['loc3']['lat']
lont = locs['20220530']['loc3']['lon']
zdr_column_height_22May30_tor3 = getBeamHeight(latt, lont, tilt[4], radar_id='CASDR')

zdr_col_hgts_22May30 = np.array([zdr_column_height_22May30_tor1, zdr_column_height_22May30_tor2, zdr_column_height_22May30_tor3])

print("---ZDR column heights, 2022-05-30---")
print(f"Tornado 1: {zdr_column_height_22May30_tor1:.3f} km")
print(f"Tornado 2: {zdr_column_height_22May30_tor2:.3f} km")
print(f"Tornado 3: {zdr_column_height_22May30_tor3:.3f} km")
# print(f"Mean:      {np.mean(zdr_col_hgts_22May30):.3f} km")
print(" ")




### May 21, 2022 - Tornadoes 3, 4 have possible weak ZDR columns
tilt = [0.0, 0.5025, 0.8925, 1.1895, 1.6780, 2.0955,
        2.7960, 3.6390, 4.4330, 5.2050, 6.7925,
        8.4075, 9.9760, 11.6375, 14.8345, 18.0315,
        21.2065, 24.3705]

latt = locs['20220521']['loc3']['lat']
lont = locs['20220521']['loc3']['lon']
zdr_column_height_22May21_tor3 = getBeamHeight(latt, lont, tilt[11], radar_id='CASKR')
# other locs: 44.1503, -79.0895 at 1330 EDT, tilt 10;   44.0199, -79.4470 at 1312 EDT, tilt 17
latt = 44.0199
lont = -79.4470
zdr_column_height_22May21_tor3_2 = getBeamHeight(latt, lont, tilt[17], radar_id='CASKR')
latt = 44.1503
lont = -79.0895
zdr_column_height_22May21_tor3_3 = getBeamHeight(latt, lont, tilt[10], radar_id='CASKR')

latt = locs['20220521']['loc4']['lat']
lont = locs['20220521']['loc4']['lon']
zdr_column_height_22May21_tor4 = getBeamHeight(latt, lont, tilt[8], radar_id='CASKR')
# other locs: 44.1640, -78.9185 at 1336 EDT, tilt 9;   44.2444, -78.6499 at 1348 EDT, tilt 7;   44.3106, -78.5093 at 1354 EDT, tilt 6
latt = 44.1640
lont = -78.9185
zdr_column_height_22May21_tor4_2 = getBeamHeight(latt, lont, tilt[9], radar_id='CASKR')
latt = 44.2444
lont = -78.6499
zdr_column_height_22May21_tor4_3 = getBeamHeight(latt, lont, tilt[7], radar_id='CASKR')
latt = 44.3106
lont = -78.5093
zdr_column_height_22May21_tor4_4 = getBeamHeight(latt, lont, tilt[6], radar_id='CASKR')

zdr_col_hgts_22May21 = np.array([zdr_column_height_22May21_tor3, zdr_column_height_22May21_tor3_2, zdr_column_height_22May21_tor3_3,
                                 zdr_column_height_22May21_tor4, zdr_column_height_22May21_tor4_2, zdr_column_height_22May21_tor4_3, zdr_column_height_22May21_tor4_4])

print("---ZDR column heights, 2022-05-21---")
print(f"Tornado 3: {zdr_column_height_22May21_tor3:.3f} km (1324 EDT)")
print(f"           {zdr_column_height_22May21_tor3_2:.3f} km (1312 EDT)")
print(f"           {zdr_column_height_22May21_tor3_3:.3f} km (1330 EDT)")
print(f"Tornado 4: {zdr_column_height_22May21_tor4:.3f} km (1342 EDT)")
print(f"           {zdr_column_height_22May21_tor4_2:.3f} km (1336 EDT)")
print(f"           {zdr_column_height_22May21_tor4_3:.3f} km (1348 EDT)")
print(f"           {zdr_column_height_22May21_tor4_4:.3f} km (1354 EDT)")
# print(f"Mean:      {np.mean(zdr_col_hgts_22May21):.3f} km")
print(" ")





### June 30, 2026 - Tornadoes 1, 2 have ZDR columns
tilt = [0.0, 0.4, 0.8, 1.2, 1.6, 2.1, 2.7, 3.5, 4.4, 5.3, 6.4, 7.7, 9.4, 11.3, 13.7, 16.6, 20.2, 24.4]

latt = locs['20260630']['loc1']['lat']
lont = locs['20260630']['loc1']['lon']
zdr_column_height_26Jun30_tor1 = getBeamHeight(latt, lont, tilt[8], radar_id='CASFT')

latt = locs['20260630']['loc2']['lat']
lont = locs['20260630']['loc2']['lon']
zdr_column_height_26Jun30_tor2 = getBeamHeight(latt, lont, tilt[7], radar_id='CASFT')

zdr_col_hgts_26Jun30 = np.array([zdr_column_height_26Jun30_tor1, zdr_column_height_26Jun30_tor2])

print("---ZDR column heights, 2026-06-30---")
print(f"Tornado 1: {zdr_column_height_26Jun30_tor1:.3f} km")
print(f"Tornado 2: {zdr_column_height_26Jun30_tor2:.3f} km")
# print(f"Mean:      {np.mean(zdr_col_hgts_26Jun30):.3f} km")
print(" ")


print("---Mean ZDR column heights---")
print(f"2021-08-11: {np.mean(zdr_col_hgts_21Aug11):.3f} km")
print(f"2025-06-23: {np.mean(zdr_col_hgts_25Jun23):.3f} km")
print(f"2022-05-30: {np.mean(zdr_col_hgts_22May30):.3f} km")
print(f"2022-05-21: {np.mean(zdr_col_hgts_22May21):.3f} km")
print(f"2026-06-30: {np.mean(zdr_col_hgts_26Jun30):.3f} km")







