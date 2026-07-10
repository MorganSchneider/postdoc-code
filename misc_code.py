# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 12:50:58 2026

@author: mschne28
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

#%% Check tendency nudging

S0 = 0.010
t1 = 0
t2 = 600
t3 = 601
t4 = 1800
i1 = int(t1)
i2 = int(t2)
i3 = int(t3)
i4 = int(t4)


time = np.linspace(0,t4,t4+1)


thp = np.zeros(shape=(len(time),), dtype=float)

thp[i1:i2] = S0/2 * (time[i1:i2]-t1)**2/(t2-t1)
thp[i2:i3] = S0 * (time[i2:i3] - t2) + thp[i2-1]
thp[i3:i4] = S0 * ((time[i3:i4]-t3)/(t4-t3)) * (t4 - (time[i3:i4]+t3)/2) + thp[i3-1]
thp[i4:] = thp[i4:] + thp[i4-1]

if S0 > 0:
    print(f"Max thpert = {np.max(thp):.2f} K at {t4} s")
elif S0 < 0:
    print(f"Min thpert = {np.min(thp):.2f} K at {t4} s")


fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')
ax.plot(time, thp, '-k', linewidth=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax.xaxis.set_major_locator(MultipleLocator(600))
ax.xaxis.set_minor_locator(MultipleLocator(60))
ax.yaxis.set_major_locator(MultipleLocator(3))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.set_xlim([0,t4])
if S0 > 0:
    ax.set_ylim([0,15])
elif S0 < 0:
    ax.set_ylim([-15,0])

plt.show()


#%% CM1 original stretched x grid

ngxy = 3

dx_inner = 125.0
dx_outer = 1875.0
nos_x_len = 100000.0
tot_x_len = 200000.0

nx = 900
dx = 500

xfl = np.zeros((4,))
xfr = np.zeros((4,))
xf = np.zeros((nx+1,))
xh = np.zeros((nx,))
uh = np.zeros((nx,))


nominal_dx = 0.5 * (dx_inner + dx_outer)

ni1 = int(np.round( (tot_x_len - nos_x_len) * 0.5 / nominal_dx ))
ni2 = int(np.round( nos_x_len / dx_inner ))
ni3 = ni1

c2 = (nominal_dx - dx_inner) / (nominal_dx**2 * (ni3-1))
c1 = (dx_inner / nominal_dx) - c2*nominal_dx

for i in range(ni1+1, ni1+ni2+3):
    tem = ni1*nominal_dx + (i-ni1-1)*dx_inner
    if (i>=1) & (i<=nx+1):
        xf[i-1] = tem
    elif i < 1:
        xfl[i-1+ngxy] = tem
    elif i > nx+1:
        xfr[i-1-(nx+1)] = tem

for i in range(ni1+ni2+2, ni1+ni2+ni3+ngxy+3):
    # tem = ni1*nominal_dx + (ni1+ni2+1-ni1-1)*dx_inner + (c1 + c2*(i-1-ni1-ni2)*nominal_dx)*(i-1-ni1-ni2)*nominal_dx
    tem = ni1*nominal_dx + ni2*dx_inner + c1*(i-1-ni1-ni2)*nominal_dx + c2*((i-1-ni1-ni2)*nominal_dx)**2
    if (i>=1) & (i<=nx+1):
        xf[i-1] = tem
    elif i < 1:
        xfl[i-1+ngxy] = tem
    elif i > nx+1:
        xfr[i-1-(nx+1)] = tem

for i in range(1-ngxy, ni1+1):
    # tem = ni1*nominal_dx + (ni1+1-ni1-1)*dx_inner - (c1 + c2*(ni1+1-i)*nominal_dx)*(ni1+1-i)*nominal_dx
    tem = ni1*nominal_dx - c1*(ni1+1-i)*nominal_dx - c2*((ni1+1-i)*nominal_dx)**2
    if (i>=1) & (i<=nx+1):
        xf[i-1] = tem
    elif i < 1:
        xfl[i-1+ngxy] = tem
    elif i > nx+1:
        xfr[i-1-(nx+1)] = tem

for i in range(1,ngxy+1):
    xfl[ngxy-i+1] = xf[0] - i*dx_outer
    xfr[i-1] = xf[nx] + i*dx_outer


for i in range(nx):
    if i == nx:
        tem = 0.5*(xfr[0] + xf[i])
    else:
        tem = 0.5*(xf[i+1] + xf[i])
    xh[i] = tem

fig,ax = plt.subplots(figsize=(8,5))
ax.scatter(xh/1000, np.zeros((nx,)), marker='.', s=1, c='k')
plt.show()


#%% Test new CM1 x grid stretching with double mesh

# known to work--
# 20,180,820 m / 10,80,180 km / nx=1400
# 25,225,2275 m / 10,100,200 km / nx=1200
# 30,270,1730 m / 12,87,187 km / nx=1000
# 30,270,2230 m / 12,87,212 km / nx=1000

ngxy = 3

dx_inner = 25.0
dx_mid = 175.0
dx_outer = 1825.0

nos_x_len = 10000.0
mid_x_len = 80000.0
tot_x_len = 180000.0

mult = 0.5

nominal_dx_inner = 0.5*(dx_inner + dx_mid)
nominal_dx_outer = 0.5*(dx_mid + dx_outer)

nx1 = nos_x_len / dx_inner
nx2 = (mid_x_len - nos_x_len) / nominal_dx_inner
nx3 = (tot_x_len - mid_x_len) / nominal_dx_outer
nx = int(nx1 + nx2 + nx3)

print(f"nx1={nx1}, nx2={nx2}, nx3={nx3}, total nx={nx1+nx2+nx3}")


ni1 = int( round( (tot_x_len-mid_x_len)*0.5/nominal_dx_outer ) )
ni2 = int( round( (mid_x_len-nos_x_len)*0.5/nominal_dx_inner ) )
ni3 = int( round( nos_x_len / dx_inner ) )
ni4 = ni2
ni5 = ni1

c2 = (nominal_dx_inner - dx_inner) / (nominal_dx_inner**2 * (ni2-1))
c1 = (dx_inner / nominal_dx_inner) - c2*nominal_dx_inner

c4 = (nominal_dx_outer - dx_mid) / (nominal_dx_outer**2 * (ni1-1))
c3 = (dx_mid / nominal_dx_outer) - c4*nominal_dx_outer


xfl = np.zeros((4,))
xfr = np.zeros((4,))
xf = np.zeros((nx+1,))
xh = np.zeros((nx,))
uh = np.zeros((nx,))



for i in range(1-ngxy, ni1+1):
    tem = ni1*nominal_dx_outer - (c3 + c4*(ni1+1-i)*nominal_dx_outer)*(ni1+1-i)*nominal_dx_outer - 0.5*tot_x_len
    if (i>=1) & (i<=nx+1):
        xf[i-1] = tem
    elif i < 1:
        xfl[i-1+ngxy] = tem
    elif i > nx+1:
        xfr[i-1-(nx+1)] = tem

for i in range(ni1+1, ni1+ni2+2):
    tem = ni1*nominal_dx_outer + ni2*nominal_dx_inner - (c1 + c2*(ni1+ni2+1-i)*nominal_dx_inner)*(ni1+ni2+1-i)*nominal_dx_inner - 0.5*tot_x_len
    if (i>=1) & (i<=nx+1):
        xf[i-1] = tem
    elif i < 1:
        xfl[i-1+ngxy] = tem
    elif i > nx+1:
        xfr[i-1-(nx+1)] = tem

for i in range(ni1+ni2+2, ni1+ni2+ni3+2):
    tem = ni1*nominal_dx_outer + ni2*nominal_dx_inner + (i-ni1-ni2-1)*dx_inner - 0.5*tot_x_len
    if (i>=1) & (i<=nx+1):
        xf[i-1] = tem
    elif i < 1:
        xfl[i-1+ngxy] = tem
    elif i > nx+1:
        xfr[i-1-(nx+1)] = tem

for i in range(ni1+ni2+ni3+2, ni1+ni2+ni3+ni4+2):
    tem = ni1*nominal_dx_outer + ni2*nominal_dx_inner + ni3*dx_inner + (c1 + c2*(i-1-ni1-ni2-ni3)*nominal_dx_inner)*(i-1-ni1-ni2-ni3)*nominal_dx_inner - 0.5*tot_x_len
    if (i>=1) & (i<=nx+1):
        xf[i-1] = tem
    elif i < 1:
        xfl[i-1+ngxy] = tem
    elif i > nx+1:
        xfr[i-1-(nx+1)] = tem

for i in range(ni1+ni2+ni3+ni4+2, ni1+ni2+ni3+ni4+ni5+ngxy+2):
    tem = ni1*nominal_dx_outer + (ni2+ni4)*nominal_dx_inner + ni3*dx_inner + (c3+c4*(i-1-ni1-ni2-ni3-ni4)*nominal_dx_outer)*(i-1-ni1-ni2-ni3-ni4)*nominal_dx_outer - 0.5*tot_x_len
    if (i>=1) & (i<=nx+1):
        xf[i-1] = tem
    elif i < 1:
        xfl[i-1+ngxy] = tem
    elif i > nx+1:
        xfr[i-1-(nx+1)] = tem


for i in range(1,ngxy+1):
    xfl[ngxy-i+1] = xf[0] - i*dx_outer
    xfr[i-1] = xf[nx] + i*dx_outer

for i in range(nx):
    if i == nx:
        tem = 0.5*(xfr[0] + xf[i])
    else:
        tem = 0.5*(xf[i+1] + xf[i])
    xh[i] = tem

yh = xh


fig,ax = plt.subplots(3, 1, figsize=(8,9))
ax[0].scatter(xf/1000, np.zeros((nx+1,)), marker='.', s=1, c='k')
ax[0].set_title('Full x grid')
ax[0].set_xlim([-tot_x_len/2/1000, tot_x_len/2/1000])

ax[1].scatter(xf/1000, np.zeros((nx+1,)), marker='.', s=1, c='k')
ax[1].set_title('Inner x mesh')
ax[1].set_xlim([0, mid_x_len/2/1000])

ax[2].scatter(xf/1000, np.zeros((nx+1,)), marker='.', s=1, c='k')
ax[2].set_title('Outer x mesh')
ax[2].set_xlim([mid_x_len/2/1000, tot_x_len/2/1000])

plt.show()



if False:
    xgrid = list(xh)
    ygrid = list(yh)
    np.savetxt("C:/Users/mschne28/Documents/input_grid_x", xgrid, fmt='%f')
    np.savetxt("C:/Users/mschne28/Documents/input_grid_y", ygrid, fmt='%f')


#%%

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import re


ds = nc.Dataset('C:/Users/mschne28/Documents/cm1out/double_stretch_test/cm1out_000001.nc')
xh = ds.variables['xh'][:].data
xf = ds.variables['xf'][:].data
ds.close()


patx = re.compile(r"stretch_x")
paty = re.compile(r"stretch_y")
with open('C:/Users/mschne28/Documents/cm1out/double_stretch_test/namelist.input', 'r') as f:
    lines = f.readlines()
l = []
for i in range(len(lines)):
    matchx = patx.search(lines[i])
    matchy = paty.search(lines[i])
    if matchx or matchy:
        l.append(lines[i+1:i+7])

strx = {}
stry = {}

underscore = re.escape("_")
pat = rf"[A-Za-z{underscore}]"

for xl in l[0]:
    varstr = "".join(re.findall(pat, xl.strip()))
    val = float(re.findall(r'\d+\.?\d+?', xl.strip())[0])
    print(f"{varstr} = {val}")
    strx.update({f"{varstr}":val})

for yl in l[1]:
    varstr = "".join(re.findall(pat, yl.strip()))
    val = float(re.findall(r'\d+\.?\d+?', yl.strip())[0])
    print(f"{varstr} = {val}")
    stry.update({f"{varstr}":val})



nx = len(xh)

i1 = np.where(np.diff(xf*1000) == 25)[0][-1]
i2 = np.where(np.diff(xf*1000) == 225)[0][-1]

fig,ax = plt.subplots(3, 1, figsize=(8,9))
ax[0].scatter(xf, np.zeros((nx+1,)), marker='.', s=1, c='k')
ax[0].set_title('Full x grid')
ax[0].set_xlim([xf[0], xf[-1]])

ax[1].scatter(xf, np.zeros((nx+1,)), marker='.', s=1, c='k')
ax[1].set_title('Inner+middle x mesh')
ax[1].set_xlim([xf[i1], xf[i2]])

ax[2].scatter(xf, np.zeros((nx+1,)), marker='.', s=1, c='k')
ax[2].set_title('Outer x mesh')
ax[2].set_xlim([xf[i2], xf[-1]])





fig,ax = plt.subplots(1, 1, figsize=(8,4))
ax.axhline(strx['dx_inner'], color='dimgray', linestyle='--', linewidth=0.75)
ax.axhline(strx['dx_mid'], color='dimgray', linestyle='--', linewidth=0.75)
ax.axhline(strx['dx_outer'], color='dimgray', linestyle='--', linewidth=0.75)
ax.axvline(x=strx['tot_x_len']/-2000, color='dimgray', linewidth=0.75, linestyle='--')
ax.axvline(x=strx['mid_x_len']/-2000, color='dimgray', linewidth=0.75, linestyle='--')
ax.axvline(x=strx['nos_x_len']/-2000, color='dimgray', linewidth=0.75, linestyle='--')
ax.axvline(x=strx['nos_x_len']/2000, color='dimgray', linewidth=0.75, linestyle='--')
ax.axvline(x=strx['mid_x_len']/2000, color='dimgray', linewidth=0.75, linestyle='--')
ax.axvline(x=strx['tot_x_len']/2000, color='dimgray', linewidth=0.75, linestyle='--')
ax.plot(xh, np.diff(xf)*1000, '-k', linewidth=2)
ax.set_xlabel('Position (km)', fontsize=12)
ax.set_ylabel('Grid spacing (m)', fontsize=12)
ax.set_ylim([0,2500])

ax.fill_between([strx['nos_x_len']/-2000, strx['nos_x_len']/2000], [0,0], [strx['dx_inner'],strx['dx_inner']],
                color='red',label="Inner mesh")
ax.fill_between([strx['mid_x_len']/-2000, strx['nos_x_len']/-2000], [strx['dx_inner'],strx['dx_inner']], [strx['dx_mid'],strx['dx_mid']],
                color='gold', alpha=0.5, label="Middle mesh")
ax.fill_between([strx['nos_x_len']/2000, strx['mid_x_len']/2000], [strx['dx_inner'],strx['dx_inner']], [strx['dx_mid'],strx['dx_mid']],
                color='gold', alpha=0.5, label="Middle mesh")
ax.fill_between([strx['tot_x_len']/-2000, strx['mid_x_len']/-2000], [strx['dx_mid'],strx['dx_mid']], [strx['dx_outer'],strx['dx_outer']],
                color='b', alpha=0.3, label="Outer mesh")
ax.fill_between([strx['mid_x_len']/2000, strx['tot_x_len']/2000], [strx['dx_mid'],strx['dx_mid']], [strx['dx_outer'],strx['dx_outer']],
                color='b', alpha=0.3, label="Outer mesh")

plt.show()



#%%

from CM1utils import *


ds = nc.Dataset('D:/cwe/semislip_wk_500m/cm1out_000001.nc')
z = ds.variables['zh'][:].data
u = ds.variables['u0'][:].data[0,:,0,0]
v = ds.variables['v0'][:].data[0,:,0,0]
th = ds.variables['th0'][:].data[0,:,0,0]
prs = ds.variables['prs0'][:].data[0,:,0,0]
qv = ds.variables['qv0'][:].data[0,:,0,0]
umove = ds.variables['umove'][:].data[0]
vmove = ds.variables['vmove'][:].data[0]
ds.close()

u = (u+umove)*units('m/s')
v = (v+vmove)*units('m/s')


T = th * (prs/100000.)**0.286
e = (qv * prs/100) / (0.622+qv)
Td = 243.5 / ((17.67/(np.log(e/6.112)))-1) + 273.15

T_parcel = mc.parcel_profile(prs/100*units.hPa, T[0]*units.K, Td[0]*units.K).to('degC')

lcl_pressure,lcl_temperature = mc.lcl(prs[0]/100*units.hPa, T[0]*units.K, Td[0]*units.K)



prs = (prs/100)*units.hPa
T = (T-273.15)*units.degC
Td = (Td-273.15)*units.degC


u_kts = u.to('kts')
v_kts = v.to('kts')




fig = plt.figure(figsize=(9,9))
skew = SkewT(fig)

skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 40)

# Plot special lines
skew.ax.grid(linewidth=1)
skew.ax.axvline(0, color='c', linestyle='--', linewidth=2) # Plot a zero degree isotherm
skew.plot_dry_adiabats(linewidth=1.5)
skew.plot_moist_adiabats(linewidth=1.5)
skew.plot_mixing_lines(linewidth=1.5)
skew.ax.axhline(y=lcl_pressure, xmin=-50, xmax=30, color='k', linestyle='--', linewidth=1.5)



skew.plot(prs, T, 'r', linewidth=2.5)
skew.plot(prs, Td, 'g', linewidth=2.5)
skew.plot(prs, T_parcel, 'k', linewidth=3) # Plot the parcel profile as a black line
skew.shade_cin(prs, T, T_parcel, alpha=0.5) # Shade areas of CAPE and CIN
skew.shade_cape(prs, T, T_parcel, alpha=0.5)
# skew.plot_barbs(prs[::6], u_kts[::6], v_kts[::6], xloc=1.1, plot_units=units('kts'))
skew.plot(lcl_pressure, lcl_temperature.to('degC'), 'ko', markerfacecolor='k', markersize=7) # Plot LCL as black dot


skew.ax.set_xlabel("Temperature (C)", fontsize=14)
skew.ax.set_ylabel("Pressure (hPa)", fontsize=14)
skew.ax.set_title("Sample Skew-T log-P + hodograph", fontsize=14)

# Create a hodograph
# ax_hod = inset_axes(skew.ax, '45%', '45%', loc=1)
# H = Hodograph(ax_hod, component_range=40.)
# H.add_grid(increment=10)
# H.plot(u*1.5, v*2, color='k', linewidth=2.5)
# ax_hod.set_xlabel('U wind (m/s)', fontsize=11)
# ax_hod.set_ylabel('V wind (m/s)', fontsize=11)

# ax_hod.quiver(u[0]*1.5, v[0]*2, u[np.argmin(abs(z-2))]*1.5, v[np.argmin(abs(z-2))]*2, color='k', linewidth=1, scale=83)
# ax_hod.plot(u[np.argmin(abs(z-2))]*1.5, v[np.argmin(abs(z-2))]*2, 'ko', markerfacecolor='w', markersize=7)
# ax_hod.quiver(u[0]*1.5, v[0]*2, u[np.argmin(abs(z-3))]*1.5, v[np.argmin(abs(z-3))]*2, color='k', linewidth=1, scale=83)
# ax_hod.plot(u[np.argmin(abs(z-3))]*1.5, v[np.argmin(abs(z-3))]*2, 'ko', markerfacecolor='w', markersize=7)
# ax_hod.quiver(u[0]*1.5, v[0]*2, u[np.argmin(abs(z-4))]*1.5, v[np.argmin(abs(z-4))]*2, color='k', linewidth=1, scale=83)
# ax_hod.plot(u[np.argmin(abs(z-4))]*1.5, v[np.argmin(abs(z-4))]*2, 'ko', markerfacecolor='w', markersize=7)



ax_hod.xaxis.set_major_locator(MultipleLocator(10))
# ax_hod.xaxis.set_minor_locator(MultipleLocator(5))
ax_hod.yaxis.set_major_locator(MultipleLocator(10))
# ax_hod.yaxis.set_minor_locator(MultipleLocator(5))



plt.show()


#%%

from era5utils import *

fp = "C:/Users/mschne28/Documents/era5/tor_outbreaks/"

# region='Alberta'
#region='SasMan'
region='GreatLakes'
# region='Canada'

yyyyt = 2026
mmt   = 4
ddt   = 15
hht   = 8

# yyyyt = 2021; mmt = 8; ddt = 11; tor hht = 20,21 (18-22)
# yyyyt = 2025; mmt = 6; ddt = 23-24; tor hht = 20,22,1 (19-23, 0-3)
# yyyyt = 2022; mmt = 5; ddt = 30-31; tor hht = 0,1,2 (22-23, 0-3)
# yyyyt = 2022; mmt = 5; ddt = 21; tor hht = 15,17 (13-18)
# yyyyt = 2025; mmt = 7; ddt = 24; bow/db hht = 21-23 (19-23)
# yyyyt = 2026; mmt = 4; ddt = 15; tor/null hht = 5,6 (2-8)

timt="%d-%s-%sT%s:00:00.000000000" %(yyyyt,str(mmt).zfill(2),str(ddt).zfill(2),str(hht).zfill(2))


lattstart = 46.4797  #Latitude of tornado/beginning of hail swath (from NTP/NHP database)
lontstart = -83.6403 #Longitude of tornado/beginning of hail swath (from NTP/NHP database)

if (yyyyt==2021) & (mmt==8) & (ddt==11):
    lattstart = [46.4797, 46.5449, 46.5724, 46.4657, 46.3870, 46.6557]
    lontstart = [-83.6403, -83.5907, -83.2889, -83.2155, -82.7493, -82.4617]
elif (yyyyt==2025) & (mmt==6) & ((ddt==23)|(ddt==24)):
    lattstart = [48.3975, 48.4647, 48.2033, 47.9149, 46.7872, 46.7757]
    lontstart = [-75.5880, -75.4900, -73.5234, -73.5056, -70.7827, -70.4072]
elif (yyyyt==2022) & (mmt==5) & ((ddt==30)|(ddt==31)):
    lattstart = [48.6148, 48.6699, 48.9114, 49.3689, 48.6757]
    lontstart = [-93.5304, -93.5219, -93.5778, -93.0754, -92.2308]
elif (yyyyt==2022) & (mmt==5) & (ddt==21):
    lattstart = [43.0179, 42.9217, 44.1058, 44.1755]
    lontstart = [-81.2216, -81.1977, -79.1458, -78.7722]
elif (yyyyt==2025) & (mmt==7) & (ddt==24):
    lattstart = [43.4700, 44.2500, 45.0000]
    lontstart = [-81.1800, -80.5000, -79.2500]
elif (yyyyt==2026) & (mmt==4) & (ddt==15):
    lattstart = [42.2718, 42.2684, 42.2393, 42.2423]
    lontstart = [-83.7524, -83.2104, -83.0576, -82.8484]




# region='Alberta'
# region='SasMan'
# region='GreatLakes'
region='Canada'

if region=='GreatLakes':    
    [lonW,lonE,latS,latN] = [-93.0,-70.0,39.5,52.0]
    n=5
    shrinkscale = 0.6
elif region=='SasMan':
    [lonW,lonE,latS,latN] = [-114.0,-88.0,46.5,60.5]
    n=5
    shrinkscale = 0.6
elif region=='Alberta':      
    [lonW,lonE,latS,latN] = [-128.0,-102.0,46.5,60.5]
    n=5
    shrinkscale = 0.6
elif region=='Canada':
    [lonW,lonE,latS,latN] = [-132.0,-58.0,35.0,75.0]
    n=10
    shrinkscale = 0.3

figfolder=fp+'figs/'
data_preslevs=fp+f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_preslevs.nc"
data_singlevs=fp+f"era5_{yyyyt}{mmt:02.0f}{ddt:02.0f}_singlevs.nc"

data = xr.open_dataset(data_preslevs)
datas = xr.open_dataset(data_singlevs)

# def PlotMap(data, datas, plotvar, pressure_level, timt, lontstart, lattstart, region="Canada", lonlatWESN=None):
### plotvar ###
### On pressure levels-
# z:  Geopotential (m2/s2)
# pv: Potential vorticity (K m2/kg/s)
# r:  Relative humidity (%)
# q:  Specific humidity (kg/kg)
# t:  Temperature (K)
# u:  U wind (m/s)
# v:  V wind (m/s)
# w:  Vertical velocity (Pa/s)
# vo: Relative vorticity (1/s)
### On single levels-
# u10:  10-m U wind (m/s)
# v10:  10-m V wind (m/s)
# d2m:  2-m dewpoint (K)
# t2m:  2-m temperature (K)
# sp:   Surface pressure (Pa)
# u100: 100-m U wind (m/s)
# v100: 100-m V wind (m/s)
# cape: CAPE (J/kg)
# cin:  CIN (J/kg)
# z:    Geopotential of surface (m2/s2)
# tcw:  Total column water (kg/m2)


dats = datas.sel(valid_time=timt)
orog = dats['z'].values/9.81 #Geopotential height of surface (m)
u10 = dats['u10'].values*1.94384449 #10-m u wind to kts
v10 = dats['v10'].values*1.94384449 #10-m v wind to kts
u100 = dats['u100'].values*1.94384449 #100-m u wind to kts
v100 = dats['v100'].values*1.94384449 #100-m v wind to kts
d2m = dats['d2m'].values-273.15 #2-m dewpoint (C)
t2m = dats['t2m'].values-273.15 #2-m temperature (C)
psfc = dats['sp'].values/100 #Surface pressure (hPa)
cape = dats['cape'].values #CAPE
cin = dats['cin'].values #CIN
mslp = psfc * (1 - (0.0065*orog)/(t2m+273.15))**(-1*(9.81*0.0289644)/(8.31447*0.0065)) # from https://absolutepressurecalculator.com/how-to-calculate-mean-sea-level-pressure.php


datt = data.sel(valid_time=timt)
prs = datt.pressure_level.values * units.hPa
lat = datt['latitude'].values
lon = datt['longitude'].values
X,Y = np.meshgrid(lon,lat)

z = datt['z'].values/9.81 #Geopotential height (m)
t = datt['t'].values-273.15 #Temperature (C)
rh = datt['r'].values*units('percent') #Relative humidity
q = datt['q'].values*units('kg/kg') #Specific humidity
u = datt['u'].values*1.94384449 #U wind to kts
v = datt['v'].values*1.94384449 #V wind to kts
omega = datt['w'].values*units('Pa/s') #omega (Pa/s)
w = mc.vertical_velocity(omega, np.moveaxis(np.tile(prs, (len(lat),len(lon),1)), -1, 0), t*units.degC, q)
pv = datt['pv'].values
vort = datt['vo'].values

lati = np.argmin(np.abs(lat-np.mean(lattstart)))
loni = np.argmin(np.abs(lon-np.mean(lontstart)))
latt = lat[lati]
lont = lon[loni]

    
#%% Plot pressure maps

colors = ['dodgerblue','lightskyblue','cyan','mediumpurple','blueviolet','mediumvioletred']
spc_wspd = ListedColormap(colors, name="spc_wspd")

plot_params = {
    'z850':    {'cm': 'LangRainbow12',  'label': "850 mb geopotential hgt (m)",  'levels':np.arange(900,1830,30)},
    'z700':    {'cm': 'LangRainbow12',  'label': "700 mb geopotential hgt (m)",  'levels':np.arange(2100,3630,30)},
    'z500':    {'cm': 'LangRainbow12',  'label': "500 mb geopotential hgt (m)",  'levels':np.arange(4500,6630,60)},
    'z300':    {'cm': 'LangRainbow12',  'label': "300 mb geopotential hgt (m)",  'levels':np.arange(7500,9930,60)},
    't850':    {'cm': 'LangRainbow12',  'label': "850 mb temperature (C)",   'levels':np.arange(-20,42,2)},
    't700':    {'cm': 'LangRainbow12',  'label': "850 mb temperature (C)",   'levels':np.arange(-30,32,2)},
    't500':    {'cm': 'LangRainbow12',  'label': "850 mb temperature (C)",   'levels':np.arange(-40,22,2)},
    'wspd':    {'cm':  spc_wspd,        'label': "Wind speed (m s$^{-1}$)",  'levels':np.arange(40,180,20)},
    'temp':    {'cm': 'HomeyerRainbow', 'label': "Temperature (C)"},
    'dewpt':   {'cm': 'YlGnBu',         'label': "Dewpoint (C)"},
    'rh':      {'cm': 'YlGn',           'label': "Relative humidity (%)"},
    'q':       {'cm': 'YlGn',           'label': "q (g kg$^{-1}$)"},
    'vort':    {'cm': 'HomeyerRainbow', 'label': "\u03B6 (s$^{-1}$)"},
    'cape':    {'cm': 'LangRainbow12',  'label': "CAPE (J kg$^{-1}$)"},
    'srh':     {'cm': 'LangRainbow12',  'label': "SRH (m$^{2}$ s$^{-2}$)"}
}



### Upper air maps ###
for presl in [300,500,700,850]:
    
    zz = datt.sel(pressure_level=presl)['z'].values/9.81
    uu = datt.sel(pressure_level=presl)['u'].values*1.94384449
    vv = datt.sel(pressure_level=presl)['v'].values*1.94384449
    wspd = np.sqrt(uu**2 + vv**2)
    tt = datt.sel(pressure_level=presl)['t'].values-273.15
    qq = datt.sel(pressure_level=presl)['q'].values
    e = (qq * presl) / (0.622+qq)
    td = 243.5 / ((17.67/(np.log(e/6.112)))-1)
    
    
    fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()})
    ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle='-')
    ax.add_feature(cfeature.STATES, linestyle='-')
    ax.add_feature(cfeature.COASTLINE)
    
    
    if presl == 300:
        cf = ax.contourf(X, Y, wspd, levels=np.linspace(60,180,7), cmap=spc_wspd, vmin=60, vmax=180, alpha=0.35)
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Wind speed [kts]', shrink=shrinkscale)
        cbar.set_label('Wind speed [kts]', fontsize=16)  
        cbar.ax.tick_params(labelsize=12)
        cw = ax.contour(X, Y, wspd, levels=np.linspace(60,160,6), colors='b', linestyles='-', linewidths=1)
        ax.clabel(cw, inline=True, fontsize=11, fmt='%d')#.1f')
    
    if presl == 500:
        cf = ax.contourf(X, Y, wspd, levels=np.linspace(40,160,7), colors=colors, vmin=40, vmax=160, alpha=0.35)
        cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Wind speed [kts]', shrink=shrinkscale)
        cbar.set_label('Wind speed [kts]', fontsize=16)  
        cbar.ax.tick_params(labelsize=12)
        cw = ax.contour(X, Y, wspd, levels=np.linspace(40,140,6), colors='b', linestyles='-', linewidths=1)
        ax.clabel(cw, inline=True, fontsize=11, fmt='%d')#.1f')
    
    if presl == 700:
        rr = np.mean(datt.sel(pressure_level=slice(700,500))['r'].values, axis=0)
        t75 = datt.sel(pressure_level=slice(700,500))['t'].values-273.15
        z75 = datt.sel(pressure_level=slice(700,500))['z'].values/9.81/1000
        lr = -1*(t75[-1,:,:]-t75[0,:,:])/(z75[-1,:,:]-z75[0,:,:])
        
        cf = ax.contourf(X, Y, rr, levels=[70,101], colors=['limegreen'], alpha=0.3)
        # cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Rel. humidity [%]', shrink=shrinkscale)
        # cbar.set_label('Rel. humidity [%]', fontsize=16)  
        # cbar.ax.tick_params(labelsize=12)
        cr = ax.contour(X, Y, rr, levels=np.arange(70,101,10), colors='g', linestyles='-', linewidths=1)
        ax.clabel(cr, inline=True, fontsize=11, fmt='%d')#.1f')
        # ct1 = ax.contour(X, Y, tt, levels=np.arange(2,10,2), colors='r', linestyles='--', linewidths=1.5)
        # ct2 = ax.contour(X, Y, tt, levels=np.arange(10,30,2), colors='r', linestyles='--', linewidths=2.5)
        # ct3 = ax.contour(X, Y, tt, levels=np.arange(-30,2,2), colors='b', linestyles='--', linewidths=1.5)
        # ax.clabel(ct1, inline=True, fontsize=11, fmt='%d')#.1f')
        # ax.clabel(ct2, inline=True, fontsize=11, fmt='%d')#.1f')
        # ax.clabel(ct3, inline=True, fontsize=11, fmt='%d')#.1f')
        cf2 = ax.contourf(X, Y, lr, levels=[8,10], colors=['r'], alpha=0.3)
        cl = ax.contour(X, Y, lr, levels=np.linspace(6,9,7), colors='r', linewidths=[1,1.5,1,1,1,1,2.5], linestyles=['--','--','-','-','-','-','-'])
        ax.clabel(cl, inline=True, fontsize=11, fmt='%d')#.1f')
    
    if presl == 850:
        cf = ax.contourf(X, Y, td, levels=np.arange(8,30,2), colors=['limegreen'], alpha=0.3)
        # cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Dewpoint [C]', shrink=shrinkscale)
        # cbar.set_label('Dewpoint [C]', fontsize=16)  
        # cbar.ax.tick_params(labelsize=12)
        cd = ax.contour(X, Y, td, levels=np.arange(8,30,2), colors='g', linestyles='-', linewidths=1)
        ct1 = ax.contour(X, Y, tt, levels=np.arange(2,20,2), colors='r', linestyles='--', linewidths=1.5)
        ct2 = ax.contour(X, Y, tt, levels=np.arange(20,40,2), colors='r', linestyles='--', linewidths=2.5)
        ct3 = ax.contour(X, Y, tt, levels=np.arange(-30,2,2), colors='b', linestyles='--', linewidths=1.5)
        ax.clabel(cd, inline=True, fontsize=11, fmt='%d')#.1f')
        ax.clabel(ct1, inline=True, fontsize=11, fmt='%d')#.1f')
        ax.clabel(ct2, inline=True, fontsize=11, fmt='%d')#.1f')
        ax.clabel(ct3, inline=True, fontsize=11, fmt='%d')#.1f')
    
    cs = ax.contour(X, Y, zz, levels=plot_params[f"z{presl}"]['levels'], colors='black', linewidths=2)
    ax.clabel(cs, inline=True, fontsize=11, fmt='%d')#.1f')

    #Wind barbs are plotted every n gridpoint
    cuv = ax.barbs(X[::n,::n], Y[::n,::n], uu[::n,::n], vv[::n,::n])
    
    
    plt.show()



### Surface maps ###


# 2-m temperature, surface pressure, 10-m wind, 100-m wind
fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()})
ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linestyle='-')
ax.add_feature(cfeature.STATES, linestyle='-')
ax.add_feature(cfeature.COASTLINE)

cf = ax.contourf(X, Y, t2m, levels=np.arange(-30,32,2), cmap='LangRainbow12', alpha=0.8)
cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Temperature [C]', shrink=shrinkscale)
cbar.set_label('Temperature [C]', fontsize=16)  
cbar.ax.tick_params(labelsize=12)
cp = ax.contour(X, Y, mslp, levels=np.arange(900,1020,2), colors='k', linewidths=1.5, linestyles='-')
ax.clabel(cp, inline=True, fontsize=11, fmt='%d')#.1f')
cuv100 = ax.barbs(X[::n,::n], Y[::n,::n], u100[::n,::n], v100[::n,::n], color='r')
cuv10 = ax.barbs(X[::n,::n], Y[::n,::n], u10[::n,::n], v10[::n,::n], color='k')

plt.show()


# 2-m dewpoint, surface pressure, 10-m wind, 100-m wind
fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()})
ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linestyle='-')
ax.add_feature(cfeature.STATES, linestyle='-')
ax.add_feature(cfeature.COASTLINE)

cf = ax.contourf(X, Y, d2m, levels=np.arange(-10,22,2), cmap='terrain_r', alpha=0.7)
cbar = plt.colorbar(cf, ax=ax, orientation='vertical', label='Dewpoint [C]', shrink=shrinkscale)
cbar.set_label('Dewpoint [C]', fontsize=16)  
cbar.ax.tick_params(labelsize=12)
cp = ax.contour(X, Y, mslp, levels=np.arange(900,1020,2), colors='k', linewidths=1.5, linestyles='-')
ax.clabel(cp, inline=True, fontsize=11, fmt='%d')#.1f')
cuv100 = ax.barbs(X[::n,::n], Y[::n,::n], u100[::n,::n], v100[::n,::n], color='r')
cuv10 = ax.barbs(X[::n,::n], Y[::n,::n], u10[::n,::n], v10[::n,::n], color='k')

plt.show()




#%%

### Composite parameter maps

zz = datt.sel(pressure_level=500)['z'].values/9.81
u3000 = InterpolateToHeightAboveGround(z, orog, u, 3000)
v3000 = InterpolateToHeightAboveGround(z, orog, v, 3000)

# CAPE, 500 mb gph
fig,ax = plt.subplots(figsize=(20,15), subplot_kw={'projection':ccrs.PlateCarree()})
ax.set_extent([lonW, lonE, latS, latN], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linestyle='-')
ax.add_feature(cfeature.STATES, linestyle='-')
ax.add_feature(cfeature.COASTLINE)

colormap = cmocean.cm.thermal.copy()
colormap.set_over("gray")
colormap.set_under(alpha=0)
vmin = 100
vmax = 3000

# norm = BoundaryNorm(np.linspace(vmin, vmax, nlevels+1), ncolors=colormap.N)
cf = ax.contourf(X, Y, cape, levels=np.arange(vmin,vmax+1,100), cmap=colormap)
cbar = plt.colorbar(cf, ax=ax, orientation='vertical', shrink=shrinkscale)
cbar.set_label('CAPE [J/kg]', fontsize=16)  
cbar.ax.tick_params(labelsize=12)

cf2 = ax.contourf(X, Y, np.ma.masked_array(cin, cape<vmin), levels=[25,100,200], colors=['w','w','w'], alpha=0.5, hatches=[None,None,'///'], edgecolor='w')
cc = ax.contour(X, Y, np.ma.masked_array(cin, cape<vmin), levels=[25,100], colors='w', linewidths=[1,1.5], linestyles='-')
ax.clabel(cc, inline=True, fontsize=11, fmt='%d')

cs = ax.contour(X, Y, zz, levels=np.arange(4500,6630,60), colors='k', linewidths=2)#'white', linewidths=1)
ax.clabel(cs, inline=True, fontsize=11, fmt='%d')

cuv = ax.barbs(X[::n,::n], Y[::n,::n], u3000[::n,::n]-u10[::n,::n], v3000[::n,::n]-v10[::n,::n], color="k")


plt.show()






































