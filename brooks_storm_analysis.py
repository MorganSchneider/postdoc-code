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










