# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 09:52:03 2026

@author: mschne28
"""

from CM1utils import *
from metpy.plots import SkewT, Hodograph
import metpy.calc as mc
from metpy.units import units
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#%% Christoph's HRDPS sounding

# HRDPS - 20 Aug 2025 - 2100 UTC - 51.166 N, -113.135 W

pres_sfc = 907.0 #mb
th_sfc = 304.5256933473404 #K
qv_sfc = 11.144332615604222 #g/kg
u_sfc = -1.8 #m/s
v_sfc = -1.2 #m/s

zh = np.array([0.0, 67.04758114500032, 557.1748619495188, 2173.760192300211, 3410.8013272430317,
      4825.37128232195, 6482.569271498218, 8497.452278389397, 9713.865092325095,
      11183.702111285256, 13061.39801286179])
th = np.array([304.5256933473404, 304.4791531562847, 304.56903313413426, 308.32229150727426, 312.6009754765594,
      317.2483035067147, 321.5026844346983, 326.3348316303594, 333.9772164835291,
      357.2305127573019,  379.40795215960156])
qv = np.array([11.144332615604222, 10.795456719467502, 10.011003056114427, 5.27711131020631 , 2.9825857222494077,
      0.6345881189949281, 0.40111795697048724, 0.16657713120401094, 0.03976440072089758,
      0.01117479828136888, 0.006594050597500024])
u = np.array([0.0, -1.8, -3.8, -1.3, 13.9,
              16.2, 17.2, 25.0, 31.2,
              31.2, 31.2])
v = np.array([0.0, -1.2, -2.3, -0.5, 5.8,
              9.8, 14.4, 22.5, 25.1,
              25.1, 25.1])

# pres = np.zeros(shape=(len(zh),), dtype=float)
# pres[0] = pres_sfc
# for i in range(len(zh)-1):
#     tmp = mc.add_height_to_pressure(pres_sfc*units.hPa, zh[i+1]*units.m)
#     pres[i+1] = tmp.magnitude
pres = np.array([907., 900., 850., 700., 600., 500., 400., 300., 250., 200., 150.])

t = th * (pres/1000.)**0.286
e = (qv/1000. * pres) / (0.622+(qv/1000.))
td = 243.5 / ((17.67/(np.log(e/6.112)))-1) + 273.15

T_parcel = mc.parcel_profile(pres*units.hPa, t[0]*units.K, td[0]*units.K)


fig = plt.figure(figsize=(8,8))

skew = SkewT(fig=fig)
skew.plot(pres, (t-273.15), '-r', linewidth=2)
skew.plot(pres, (td-273.15), '-g', linewidth=2)
skew.plot(pres, np.array(T_parcel.magnitude[:])-273.15, '-k', linewidth=2)
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 30)
plt.title('HRDPS 20250820_21z 55.116,-113.135')
ax_hod = inset_axes(skew.ax, '42%', '42%', loc=1)
H = Hodograph(ax_hod, component_range=40.)
H.add_grid(increment=10)
H.plot(u, v, color='k', linewidth=1.5)

plt.show()




# Calculate sounding parameters
bwnd = mc.bunkers_storm_motion(pres*units.hPa, u*units('m/s'), v*units('m/s'), zh*units.m)
uBR = bwnd[0].magnitude[0]
vBR = bwnd[0].magnitude[1]
uBL = bwnd[1].magnitude[0]
vBL = bwnd[1].magnitude[1]
u06 = bwnd[2].magnitude[0]
v06 = bwnd[2].magnitude[1]
smBR = np.sqrt(uBR**2 + vBR**2)
angBR = 180 + np.arctan2(uBR, vBR)*180/np.pi
VH06 = np.sqrt(u06**2 + v06**2)
ang06 = 180 + np.arctan2(u06, v06)*180/np.pi


CC = mc.cape_cin(pres*units.hPa, t*units.K, td*units.K, T_parcel)
cape = CC[0].magnitude
cin = CC[1].magnitude
# MU = mc.most_unstable_cape_cin(pres*units.hPa, t*units.K, td*units.K, height=zh*units.m)
# mucape = MU[0].magnitude
# mucin = MU[1].magnitude
# SB = mc.surface_based_cape_cin(pres*units.hPa, t*units.K, td*units.K)
# sbcape = SB[0].magnitude
# sbcin = SB[1].magnitude
# ML = mc.mixed_layer_cape_cin(pres*units.hPa, t*units.K, td*units.K, height=zh*units.m)
# mlcape = ML[0].magnitude
# mlcin = ML[1].magnitude
# DC = mc.downdraft_cape(pres*units.hPa, t*units.K, td*units.K)
# dcape = DC[0].magnitude

# lcl = mc.lcl(pres*units.hPa, t[0]*units.K, td[0]*units.K)
# plcl = lcl[0].magnitude[0]
# ilcl = np.where(pres <= plcl)[0][0]
# zlcl = zh[ilcl]

# el = mc.el(pres*units.hPa, t*units.K, td*units.K)
# pel = el[0].magnitude
# iel = np.where(pres <= pel)[0][0]
# zel = zh[iel]


# shr6km = mc.bulk_shear(pres*units.hPa, u*units('m/s'), v*units('m/s'), height=zh*units.m, bottom=10*units.m, depth=6000*units.m)
# ushr06 = shr6km[0].magnitude
# vshr06 = shr6km[1].magnitude
# shr06 = np.sqrt(ushr06**2 + vshr06**2)
# angSHR = 180 + np.arctan2(ushr06, vshr06)*180/np.pi

# srh5 = mc.storm_relative_helicity(zh*units.m, u*units('m/s'), v*units('m/s'), depth=500*units.m, storm_u=uBR*units('m/s'), storm_v=vBR*units('m/s'))
# srh500 = srh5[2].magnitude
# srh1 = mc.storm_relative_helicity(zh*units.m, u*units('m/s'), v*units('m/s'), depth=1000*units.m, storm_u=uBR*units('m/s'), storm_v=vBR*units('m/s'))
# srh1km = srh1[2].magnitude
# srh3 = mc.storm_relative_helicity(zh*units.m, u*units('m/s'), v*units('m/s'), depth=3000*units.m, storm_u=uBR*units('m/s'), storm_v=vBR*units('m/s'))
# srh3km = srh3[2].magnitude


print("---HRDPS profile---")
print(f"Profile top:          {max(zh):.0f} m")
print(f"Bunkers RM:           {smBR:.1f} m/s at {angBR:.0f} deg (Vector: {uBR:.1f} m/s, {vBR:.1f} m/s)")
print(f"0-6 km mean wind:     {VH06:.1f} m/s at {ang06:.0f} deg (Vector: {u06:.1f} m/s, {v06:.1f} m/s)")
print(f"CAPE,CIN:             {cape:.0f} J/kg, {cin:.0f} J/kg")



#%% Correct Jack's ERA5 soundings from ASL to AGL and resave

import pandas as pd


fp = 'C:/Users/mschne28/Documents/brooks_storm/'
save_flag = False

# Most favourable profile --- 1121 m?
fn1 = 'ERA5_Profile_20250820_21Z_50.75_114.0_CM1.txt'

with open(fp+fn1, 'r') as f:
    lines = f.readlines()

prs0 = float(lines[0][0:10].strip())
th0 = float(lines[0][10:21].strip())
qv0 = float(lines[0][21:32].strip())
z = np.zeros(shape=(len(lines)-1,), dtype=float)
th = np.zeros(shape=(len(lines)-1,), dtype=float)
qv = np.zeros(shape=(len(lines)-1,), dtype=float)
u = np.zeros(shape=(len(lines)-1,), dtype=float)
v = np.zeros(shape=(len(lines)-1,), dtype=float)

z_offset = float(lines[1][0:10].strip()) - 10

for i in np.arange(1,len(lines)):
    z_asl = float(lines[i][0:10].strip())
    z[i-1] = z_asl - z_offset
    th[i-1] = float(lines[i][10:21].strip())
    qv[i-1] = float(lines[i][21:32].strip())
    u[i-1] = float(lines[i][32:43].strip())
    v[i-1] = float(lines[i][43:54].strip())

if save_flag:
    fsave = 'input_sounding_era5_50.75_114.0'
    hd = np.zeros((1,3))
    hd[0][0] = prs0
    hd[0][1] = th0
    hd[0][2] = qv0
    np.savetxt(fp+fsave, hd, fmt='%f')
    
    zsave = list(z)
    thsave = list(th)
    qvsave = list(qv)
    usave = list(u)
    vsave = list(v)
    
    dat1 = {'z':zsave, 'theta':thsave, 'qv':qvsave, 'u':usave, 'v':vsave}
    df1 = pd.DataFrame(data=dat1, dtype=float)
    with open(fp+fsave, 'a') as ff1:
        ff1.write(df1.to_string(header=False, index=False))

del lines




# Good thermo but probably won't be the right storm mode --- 1104 m ?
fn2 = 'ERA5_Profile_20250820_21Z_51.0_114.25_CM1.txt'

with open(fp+fn2, 'r') as f:
    lines = f.readlines()

prs0 = float(lines[0][0:10].strip())
th0 = float(lines[0][10:21].strip())
qv0 = float(lines[0][21:32].strip())
z = np.zeros(shape=(len(lines)-1,), dtype=float)
th = np.zeros(shape=(len(lines)-1,), dtype=float)
qv = np.zeros(shape=(len(lines)-1,), dtype=float)
u = np.zeros(shape=(len(lines)-1,), dtype=float)
v = np.zeros(shape=(len(lines)-1,), dtype=float)

z_offset = float(lines[1][0:10].strip()) - 10

for i in np.arange(1,len(lines)):
    z_asl = float(lines[i][0:10].strip())
    z[i-1] = z_asl - z_offset
    th[i-1] = float(lines[i][10:21].strip())
    qv[i-1] = float(lines[i][21:32].strip())
    u[i-1] = float(lines[i][32:43].strip())
    v[i-1] = float(lines[i][43:54].strip())

if save_flag:
    fsave = 'input_sounding_era5_51.0_114.25'
    hd = np.zeros((1,3))
    hd[0][0] = prs0
    hd[0][1] = th0
    hd[0][2] = qv0
    np.savetxt(fp+fsave, hd, fmt='%f')
    
    zsave = list(z)
    thsave = list(th)
    qvsave = list(qv)
    usave = list(u)
    vsave = list(v)
    
    dat2 = {'z':zsave, 'theta':thsave, 'qv':qvsave, 'u':usave, 'v':vsave}
    df2 = pd.DataFrame(data=dat2, dtype=float)
    with open(fp+fsave, 'a') as ff2:
        ff2.write(df2.to_string(header=False, index=False))

prs02 = prs0
th02 = th0
qv02 = qv0
del lines




# Closest to Christoph's HRDPS profile --- 1034 m ?
fn3 = 'ERA5_Profile_20250820_21Z_51.25_113.25_CM1.txt'

with open(fp+fn3, 'r') as f:
    lines = f.readlines()

prs0 = float(lines[0][0:10].strip())
th0 = float(lines[0][10:21].strip())
qv0 = float(lines[0][21:32].strip())
z = np.zeros(shape=(len(lines)-1,), dtype=float)
th = np.zeros(shape=(len(lines)-1,), dtype=float)
qv = np.zeros(shape=(len(lines)-1,), dtype=float)
u = np.zeros(shape=(len(lines)-1,), dtype=float)
v = np.zeros(shape=(len(lines)-1,), dtype=float)

z_offset = float(lines[1][0:10].strip()) - 10

for i in np.arange(1,len(lines)):
    z_asl = float(lines[i][0:10].strip())
    z[i-1] = z_asl - z_offset
    th[i-1] = float(lines[i][10:21].strip())
    qv[i-1] = float(lines[i][21:32].strip())
    u[i-1] = float(lines[i][32:43].strip())
    v[i-1] = float(lines[i][43:54].strip())

if save_flag:
    fsave = 'input_sounding_era5_51.25_113.25'
    hd = np.zeros((1,3))
    hd[0][0] = prs0
    hd[0][1] = th0
    hd[0][2] = qv0
    np.savetxt(fp+fsave, hd, fmt='%f')
    
    zsave = list(z)
    thsave = list(th)
    qvsave = list(qv)
    usave = list(u)
    vsave = list(v)
    
    dat3 = {'z':zsave, 'theta':thsave, 'qv':qvsave, 'u':usave, 'v':vsave}
    df3 = pd.DataFrame(data=dat3, dtype=float)
    with open(fp+fsave, 'a') as ff3:
        ff3.write(df3.to_string(header=False, index=False))

prs03 = prs0
th03 = th0
qv03 = qv0
del lines



#%% Get ERA5 sounding parameters


fp = 'C:/Users/mschne28/Documents/brooks_storm/'
fn1 = 'input_sounding_era5_50.75_114.0'
fn2 = 'input_sounding_era5_51.0_114.25'
fn3 = 'input_sounding_era5_51.25_113.25'


with open(fp+fn1, 'r') as f:
    lines = f.readlines()
prs01 = float(lines[0][0:10].strip())
th01 = float(lines[0][10:21].strip())
qv01 = float(lines[0][21:31].strip())
z1 = np.zeros(shape=(len(lines)-1,), dtype=float)
th1 = np.zeros(shape=(len(lines)-1,), dtype=float)
qv1 = np.zeros(shape=(len(lines)-1,), dtype=float)
u1 = np.zeros(shape=(len(lines)-1,), dtype=float)
v1 = np.zeros(shape=(len(lines)-1,), dtype=float)
for i in np.arange(1,len(lines)):
    z1[i-1] = float(lines[i][0:8].strip())
    th1[i-1] = float(lines[i][8:16].strip())
    qv1[i-1] = float(lines[i][16:21].strip())
    u1[i-1] = float(lines[i][21:30].strip())
    v1[i-1] = float(lines[i][30:38].strip())


with open(fp+fn2, 'r') as f:
    lines = f.readlines()
prs02 = float(lines[0][0:10].strip())
th02 = float(lines[0][10:21].strip())
qv02 = float(lines[0][21:31].strip())
z2 = np.zeros(shape=(len(lines)-1,), dtype=float)
th2 = np.zeros(shape=(len(lines)-1,), dtype=float)
qv2 = np.zeros(shape=(len(lines)-1,), dtype=float)
u2 = np.zeros(shape=(len(lines)-1,), dtype=float)
v2 = np.zeros(shape=(len(lines)-1,), dtype=float)
for i in np.arange(1,len(lines)):
    z2[i-1] = float(lines[i][0:8].strip())
    th2[i-1] = float(lines[i][8:16].strip())
    qv2[i-1] = float(lines[i][16:22].strip())
    u2[i-1] = float(lines[i][22:31].strip())
    v2[i-1] = float(lines[i][31:39].strip())


with open(fp+fn3, 'r') as f:
    lines = f.readlines()
prs03 = float(lines[0][0:10].strip())
th03 = float(lines[0][10:21].strip())
qv03 = float(lines[0][21:31].strip())
z3 = np.zeros(shape=(len(lines)-1,), dtype=float)
th3 = np.zeros(shape=(len(lines)-1,), dtype=float)
qv3 = np.zeros(shape=(len(lines)-1,), dtype=float)
u3 = np.zeros(shape=(len(lines)-1,), dtype=float)
v3 = np.zeros(shape=(len(lines)-1,), dtype=float)
for i in np.arange(1,len(lines)):
    z3[i-1] = float(lines[i][0:8].strip())
    th3[i-1] = float(lines[i][8:16].strip())
    qv3[i-1] = float(lines[i][16:21].strip())
    u3[i-1] = float(lines[i][21:30].strip())
    v3[i-1] = float(lines[i][30:38].strip())






prs1 = np.zeros(shape=(len(z1),), dtype=float)
# prs1[0] = prs01
for i in range(len(z1)):
    tmp = mc.add_height_to_pressure(prs01*units.hPa, z1[i]*units.m)
    prs1[i] = tmp.magnitude

prs2 = np.zeros(shape=(len(z2),), dtype=float)
# prs2[0] = prs02
for i in range(len(z2)):
    tmp = mc.add_height_to_pressure(prs02*units.hPa, z2[i]*units.m)
    prs2[i] = tmp.magnitude

prs3 = np.zeros(shape=(len(z3),), dtype=float)
# prs3[0] = prs03
for i in range(len(z3)):
    tmp = mc.add_height_to_pressure(prs03*units.hPa, z3[i]*units.m)
    prs3[i] = tmp.magnitude


i1 = np.where(prs1 == np.nanmin(prs1))[0][0]
i2 = np.where(prs2 == np.nanmin(prs2))[0][0]
i3 = np.where(prs3 == np.nanmin(prs3))[0][0]

prs1 = prs1[:i1]
z1 = z1[:i1]
th1 = th1[:i1]
qv1 = qv1[:i1]
u1 = u1[:i1]
v1 = v1[:i1]

prs2 = prs2[:i2]
z2 = z2[:i2]
th2 = th2[:i2]
qv2 = qv2[:i2]
u2 = u2[:i2]
v2 = v2[:i2]

prs3 = prs3[:i3]
z3 = z3[:i3]
th3 = th3[:i3]
qv3 = qv3[:i3]
u3 = u3[:i3]
v3 = v3[:i3]



T1 = th1 * (prs1/1000.)**0.286
e = (qv1/1000. * prs1) / (0.622+(qv1/1000.))
Td1 = 243.5 / ((17.67/(np.log(e/6.112)))-1) + 273.15
T1_parcel = mc.parcel_profile(prs1*units.hPa, T1[0]*units.K, Td1[0]*units.K)

T2 = th2 * (prs2/1000.)**0.286
e = (qv2/1000. * prs2) / (0.622+(qv2/1000.))
Td2 = 243.5 / ((17.67/(np.log(e/6.112)))-1) + 273.15
T2_parcel = mc.parcel_profile(prs2*units.hPa, T2[0]*units.K, Td2[0]*units.K)

T3 = th3 * (prs3/1000.)**0.286
e = (qv3/1000. * prs3) / (0.622+(qv3/1000.))
Td3 = 243.5 / ((17.67/(np.log(e/6.112)))-1) + 273.15
T3_parcel = mc.parcel_profile(prs3*units.hPa, T3[0]*units.K, Td3[0]*units.K)






# Calculate sounding parameters
bwnd = mc.bunkers_storm_motion(prs1*units.hPa, u1*units('m/s'), v1*units('m/s'), z1*units.m)
uBR_1 = bwnd[0].magnitude[0]
vBR_1 = bwnd[0].magnitude[1]
u06_1 = bwnd[2].magnitude[0]
v06_1 = bwnd[2].magnitude[1]
smBR_1 = np.sqrt(uBR_1**2 + vBR_1**2)
angBR_1 = 180 + np.arctan2(uBR_1, vBR_1)*180/np.pi
VH06_1 = np.sqrt(u06_1**2 + v06_1**2)
ang06_1 = 180 + np.arctan2(u06_1, v06_1)*180/np.pi
CC = mc.cape_cin(prs1*units.hPa, T1*units.K, Td1*units.K, T1_parcel)
cape1 = CC[0].magnitude
cin1 = CC[1].magnitude

bwnd = mc.bunkers_storm_motion(prs2*units.hPa, u2*units('m/s'), v2*units('m/s'), z2*units.m)
uBR_2 = bwnd[0].magnitude[0]
vBR_2 = bwnd[0].magnitude[1]
u06_2 = bwnd[2].magnitude[0]
v06_2 = bwnd[2].magnitude[1]
smBR_2 = np.sqrt(uBR_2**2 + vBR_2**2)
angBR_2 = 180 + np.arctan2(uBR_2, vBR_2)*180/np.pi
VH06_2 = np.sqrt(u06_2**2 + v06_2**2)
ang06_2 = 180 + np.arctan2(u06_2, v06_2)*180/np.pi
CC = mc.cape_cin(prs2*units.hPa, T2*units.K, Td2*units.K, T2_parcel)
cape2 = CC[0].magnitude
cin2 = CC[1].magnitude

bwnd = mc.bunkers_storm_motion(prs3*units.hPa, u3*units('m/s'), v3*units('m/s'), z3*units.m)
uBR_3 = bwnd[0].magnitude[0]
vBR_3 = bwnd[0].magnitude[1]
u06_3 = bwnd[2].magnitude[0]
v06_3 = bwnd[2].magnitude[1]
smBR_3 = np.sqrt(uBR_3**2 + vBR_3**2)
angBR_3 = 180 + np.arctan2(uBR_3, vBR_3)*180/np.pi
VH06_3 = np.sqrt(u06_3**2 + v06_3**2)
ang06_3 = 180 + np.arctan2(u06_3, v06_3)*180/np.pi
CC = mc.cape_cin(prs3*units.hPa, T3*units.K, Td3*units.K, T3_parcel)
cape3 = CC[0].magnitude
cin3 = CC[1].magnitude


print("---ERA5 profiles---")
print(f"   50.75, -114.0")
print(f"Bunkers RM:           {smBR_1:.1f} m/s at {angBR_1:.0f} deg (Vector: {uBR_1:.1f} m/s, {vBR_1:.1f} m/s)")
print(f"0-6 km mean wind:     {VH06_1:.1f} m/s at {ang06_1:.0f} deg (Vector: {u06_1:.1f} m/s, {v06_1:.1f} m/s)")
print(f"CAPE,CIN:             {cape1:.0f} J/kg, {cin1:.0f} J/kg")
print(" ")
print(f"   51.0, -114.25")
print(f"Bunkers RM:           {smBR_2:.1f} m/s at {angBR_2:.0f} deg (Vector: {uBR_2:.1f} m/s, {vBR_2:.1f} m/s)")
print(f"0-6 km mean wind:     {VH06_2:.1f} m/s at {ang06_2:.0f} deg (Vector: {u06_2:.1f} m/s, {v06_2:.1f} m/s)")
print(f"CAPE,CIN:             {cape2:.0f} J/kg, {cin2:.0f} J/kg")
print(" ")
print(f"   51.25, -113.25")
print(f"Bunkers RM:           {smBR_3:.1f} m/s at {angBR_3:.0f} deg (Vector: {uBR_3:.1f} m/s, {vBR_3:.1f} m/s)")
print(f"0-6 km mean wind:     {VH06_3:.1f} m/s at {ang06_3:.0f} deg (Vector: {u06_3:.1f} m/s, {v06_3:.1f} m/s)")
print(f"CAPE,CIN:             {cape3:.0f} J/kg, {cin3:.0f} J/kg")












