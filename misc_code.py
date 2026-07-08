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


#%%




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


