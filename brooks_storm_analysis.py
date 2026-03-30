# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 15:05:24 2026

@author: mschne28
"""

from CM1utils import *

#%% Plot simulations


fp = 'C:/Users/mschne28/Documents/cm1out/brooks/era5_125m_test2/'

ds = nc.Dataset(fp+'cm1out_000029.nc')
time = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['zh'][:].data
dbz = ds.variables['dbz'][:].data[0,0,:,:]
thpert = ds.variables['th'][:].data[0,0,:,:] - ds.variables['th0'][:].data[0,0,:,:]
winterp = ds.variables['winterp'][:].data[0,:,:,:]
uinterp = ds.variables['uinterp'][:].data[0,0,:,:] + ds.variables['umove'][:].data
vinterp = ds.variables['vinterp'][:].data[0,0,:,:] + ds.variables['vmove'][:].data
zvort = ds.variables['zvort'][:].data[0,:,:,:]
umove = ds.variables['umove'][:].data[0]
vmove = ds.variables['vmove'][:].data[0]
ds.close()

iz = np.where(zh>=3)[0][0]
wmax = np.max(winterp[0:iz,:,:], axis=0)

iz = np.where(zh>=1)[0][0]
vortmax = np.max(zvort[0:iz,:,:], axis=0)



fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, dbz, 'dbz', ax, levels=np.linspace(0,70,15), datalims=[0,70],
              cmap='HomeyerRainbow', cbfs=12)
ax.contour(xh, yh, thpert, levels=[-2], colors='dodgerblue', linewidths=[1], linestyles='--')
ax.contour(xh, yh, wmax, levels=[5,10], colors=['k','k'], linewidths=[0.5,1.25])
ax.contour(xh, yh, vortmax, levels=[0.03,0.05], colors=['maroon','r'], linewidths=[1,1.5], linestyles='-')
ax.set_xlim([-100,100])
ax.set_ylim([-100,100])
ax.set_title(f"Surface dBZ (t = {time/60:.0f} min)", fontsize=12)
ax.set_xlabel("x (km)", fontsize=12)
ax.set_ylabel("y (km)", fontsize=12)

l7,= ax.plot([-150,-149], [-150,-149], color='dodgerblue', linestyle='--', linewidth=1)
# l8,= ax.plot([-150,-149], [-150,-149], color='dodgerblue', linestyle='--', linewidth=1.5)
# l1,= ax.plot([-150,-149], [-150,-149], color='#516572', linestyle='-', linewidth=1)
l1,= ax.plot([-150,-149], [-150,-149], color='k', linestyle='-', linewidth=0.5)
l2,= ax.plot([-150,-149], [-150,-149], color='k', linestyle='-', linewidth=1.25)
# l3,= ax.plot([-150,-149], [-150,-149], color='k', linestyle='-', linewidth=2.5)
l4,= ax.plot([-150,-149], [-150,-149], color='maroon', linestyle='-', linewidth=1)
l5,= ax.plot([-150,-149], [-150,-149], color='r', linestyle='-', linewidth=1.5)
# l6,= ax.plot([-150,-149], [-150,-149], color='maroon', linestyle='-', linewidth=1.5)

ax.plot([-50,-50], [-50,50], 'dimgray', [50,50], [-50,50], 'dimgray', [-50,50], [-50,-50], 'dimgray', [-50,50], [50,50], 'dimgray', linewidth=1, linestyle='--')

ax.legend(handles=[l1,l2,l4,l5,l7],
          labels=["w=5 m s$^{-1}$", "w=10 m s$^{-1}$",
                  "\u03B6=0.03 s$^{-1}$", "\u03B6=0.05 s$^{-1}$",
                  "\u03B8'=-2 K"],
          ncols=3, fontsize=10, loc='lower right')
# plt.savefig(fp+f"dbz_{time/60:.0f}min_125m_run.png", dpi=300)

plt.show()
#%%


fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, wmax, 'w', ax, levels=np.linspace(-15,15,31), datalims=[-15,15],
              cmap='balance', cbticks=np.linspace(-15,15,11), cbfs=12)
ax.contour(xh, yh, vortmax, levels=[0.025,0.05], colors=['b','dodgerblue'], linewidths=[1,1.5])
ax.set_xlim([-100,100])
ax.set_ylim([-100,100])
ax.set_title(f"Max 0-3 km w, max 0-1 km \u03B6 (t = {time/60:.0f} min)", fontsize=12)
ax.set_xlabel('x (km)', fontsize=12)
ax.set_ylabel('y (km)', fontsize=12)

l1,= ax.plot([-150,-149], [-150,-149], color='b', linestyle='-', linewidth=1)
l2,= ax.plot([-150,-149], [-150,-149], color='dodgerblue', linestyle='-', linewidth=1.5)
# l3,= ax.plot([-150,-149], [-150,-149], color='k', linestyle='-', linewidth=1)
ax.legend(handles=[l1,l2], labels=["\u03B6=0.025 s$^{-1}$", "\u03B6=0.05 s$^{-1}$"],
          fontsize=10, loc='lower right')
# plt.savefig(fp+f"wmax_{time/60:.0f}min_125m_run.png", dpi=300)



qix = 40
fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, thpert, 'thpert', ax, levels=np.linspace(-10,10,21), datalims=[-10,10],
              cmap='balance', cbticks=np.linspace(-10,10,11), cbfs=12)
ax.quiver(xh[::qix], yh[::qix], uinterp[::qix,::qix], vinterp[::qix,::qix], color='k', scale=200, width=0.003, pivot='middle')
ax.set_xlim([-100,100])
ax.set_ylim([-100,100])
ax.set_title(f"Surface \u03B8', surface winds (t = {time/60:.0f} min)", fontsize=12)
ax.set_xlabel('x (km)', fontsize=12)
ax.set_ylabel('y (km)', fontsize=12)
# plt.savefig(fp+f"thpert_{time/60:.0f}min_125m_run.png", dpi=300)



plt.show()

i_wmax = np.where(wmax == np.max(wmax))
i_vmax = np.where(vortmax == np.max(vortmax))

x_mc = np.mean([xh[i_wmax[1][0]], xh[i_vmax[1][0]]])
y_mc = np.mean([yh[i_wmax[0][0]], yh[i_vmax[0][0]]])

u_trans = x_mc*1000/time
v_trans = y_mc*1000/time

umove_new = u_trans + umove
vmove_new = v_trans + vmove



#%% Check stats time series

fp = 'C:/Users/mschne28/Documents/cm1out/brooks/era5_125m_test1/'

ds = nc.Dataset(fp+'cm1out_stats.nc')
time = ds.variables['time'][:].data
vort1km = ds.variables['vort1km'][:].data
vort3km = ds.variables['vort3km'][:].data
vort5km = ds.variables['vort5km'][:].data
wmax1000 = ds.variables['wmax1000'][:].data
wmax2500 = ds.variables['wmax2500'][:].data
wmax5000 = ds.variables['wmax5000'][:].data
ds.close()


fig,ax = plt.subplots(1, 1, figsize=(10,4), layout='constrained')
l1, = ax.plot(time/60, movmean(vort1km,5), '-r', linewidth=2)
l2, = ax.plot(time/60, movmean(vort3km,5), '-b', linewidth=2)
l3, = ax.plot(time/60, movmean(vort5km,5), '-k', linewidth=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax.xaxis.set_major_locator(MultipleLocator(30))
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(0.04))
ax.yaxis.set_minor_locator(MultipleLocator(0.01))
ax.set_xlim([0,420])
ax.set_ylim([0,0.2])
ax.legend(handles=[l1,l2,l3], labels=["\u03B6(1 km)", "\u03B6(3 km)", "\u03B6(5 km)"], loc='lower right', fontsize=12)
ax.set_title("Domain-maximum \u03B6", fontsize=14)
ax.set_xlabel('Time (min)', fontsize=12)
ax.set_ylabel("\u03B6 (1/s)", fontsize=12)
# plt.savefig(fp+'vortmax_timeseries.png', dpi=300)


fig,ax = plt.subplots(1, 1, figsize=(10,4), layout='constrained')
l1, = ax.plot(time/60, wmax1000, '-r', linewidth=2)
l2, = ax.plot(time/60, wmax2500, '-b', linewidth=2)
l3, = ax.plot(time/60, wmax5000, '-k', linewidth=2)
ax.grid(visible=True, which='major', color='dimgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax.xaxis.set_major_locator(MultipleLocator(30))
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.set_xlim([0,420])
ax.set_ylim([0,30])
ax.legend(handles=[l1,l2,l3], labels=['w(1 km)', 'w(2.5 km)', 'w(5 km)'], loc='lower right', fontsize=12)
ax.set_title('Domain-maximum w', fontsize=14)
ax.set_xlabel('Time (min)', fontsize=12)
ax.set_ylabel('w (m/s)', fontsize=12)
# plt.savefig(fp+'wmax_timeseries.png', dpi=300)

plt.show()





#%% Translated swaths

fp = 'C:/Users/mschne28/Documents/cm1out/brooks/era5_125m_test1/'

# fn = np.linspace(5,29,7)




ds = nc.Dataset(fp+'cm1out_000005.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
umove = ds.variables['umove'][:].data[0]
vmove = ds.variables['vmove'][:].data[0]
shs1 = ds.variables['shs2'][:].data[0,:,:]
hail1 = ds.variables['hail2'][:].data[0,:,:]
dbz1 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000009.nc')
shs2 = ds.variables['shs2'][:].data[0,:,:]
hail2 = ds.variables['hail2'][:].data[0,:,:]
dbz2 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000013.nc')
shs3 = ds.variables['shs2'][:].data[0,:,:]
hail3 = ds.variables['hail2'][:].data[0,:,:]
dbz3 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000017.nc')
shs4 = ds.variables['shs2'][:].data[0,:,:]
hail4 = ds.variables['hail2'][:].data[0,:,:]
dbz4 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000021.nc')
shs5 = ds.variables['shs2'][:].data[0,:,:]
hail5 = ds.variables['hail2'][:].data[0,:,:]
dbz5 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000025.nc')
shs6 = ds.variables['shs2'][:].data[0,:,:]
hail6 = ds.variables['hail2'][:].data[0,:,:]
dbz6 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000029.nc')
shs7 = ds.variables['shs2'][:].data[0,:,:]
hail7 = ds.variables['hail2'][:].data[0,:,:]
dbz7 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()


x_added = umove*3600/1000
y_added = vmove*3600/1000

xh1 = xh + 100
xh2 = xh1 + x_added
xh3 = xh1 + 2*x_added
xh4 = xh1 + 3*x_added
xh5 = xh1 + 4*x_added
xh6 = xh1 + 5*x_added
xh7 = xh1 + 6*x_added

yh = yh + 100
yh1 = yh
yh2 = yh1 + y_added
yh3 = yh1 + 2*y_added
yh4 = yh1 + 3*y_added
y_added = 2*3600/1000
yh5 = yh1 + 4*y_added
yh6 = yh1 + 5*y_added
yh7 = yh1 + 6*y_added

#%%
levs = [300,500]
cols = ['dimgray','k']
lws = [0.6,1]
dbz_levs = np.linspace(0,70,15)
dbz_cols = ['limegreen','gold','darkorange','r']
dbz_lws = [1.25,1.25,1,1]

fig,ax = plt.subplots(1, 1, figsize=(8.5,2.5), subplot_kw=dict(aspect=1), layout='constrained')

c = ax.contourf(xh1, yh, np.ma.masked_array(dbz1, dbz1<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')
ax.contourf(xh2, yh, np.ma.masked_array(dbz2, dbz2<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')
ax.contourf(xh3, yh, np.ma.masked_array(dbz3, dbz3<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')
ax.contourf(xh4, yh, np.ma.masked_array(dbz4, dbz4<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')
ax.contourf(xh5, yh, np.ma.masked_array(dbz5, dbz5<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')
ax.contourf(xh6, yh, np.ma.masked_array(dbz6, dbz6<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')
ax.contourf(xh7, yh, np.ma.masked_array(dbz7, dbz7<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')

cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_ticks(np.linspace(0,70,8))
cb.set_label('Radar reflectivity (dBZ)', fontsize=10)

ax.contour(xh1, yh, shs1, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh2, yh, shs2, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh3, yh, shs3, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh4, yh, shs4, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh5, yh, shs5, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh6, yh, shs6, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh7, yh, shs7, levels=levs, colors=cols, linewidths=lws)
ax.set_xlim([0,450])
ax.set_ylim([60,180])
ax.set_xlabel('Translated x (km)', fontsize=10)
ax.set_ylabel('y (km)', fontsize=10)
ax.set_title('Brooks test run #1: Translated updraft helicity swaths (0-7 h)', fontsize=12)
l1, = ax.plot([-2,-1], [-2,-1], 'gray', linewidth=0.75)
l2, = ax.plot([-2,-1], [-2,-1], 'k', linewidth=1)
ax.legend(handles=[l1,l2], labels=['300 m2/s2','500 m2/s2'], loc='lower right', fontsize=10)

ax.text(77, 73, '1 h', fontsize=10, fontweight='bold')
ax.text(120, 78, '2 h', fontsize=10, fontweight='bold')
ax.text(165, 80, '3 h', fontsize=10, fontweight='bold')
ax.text(215, 83, '4 h', fontsize=10, fontweight='bold')
ax.text(265, 89, '5 h', fontsize=10, fontweight='bold')
ax.text(320, 101, '6 h', fontsize=10, fontweight='bold')
ax.text(380, 102, '7 h', fontsize=10, fontweight='bold')

# plt.savefig(fp+'dbz_uh_swath_translated.png', dpi=300)

#%%

levs = [0.01,0.1]
cols = ['gray','k']
lws = [1,1]

fig,ax = plt.subplots(1, 1, figsize=(8.5,2.5), subplot_kw=dict(aspect=1), layout='constrained')

c = ax.contourf(xh1, yh, np.ma.masked_array(dbz1, dbz1<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')
ax.contourf(xh2, yh, np.ma.masked_array(dbz2, dbz2<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')
ax.contourf(xh3, yh, np.ma.masked_array(dbz3, dbz3<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')
ax.contourf(xh4, yh, np.ma.masked_array(dbz4, dbz4<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')
ax.contourf(xh5, yh, np.ma.masked_array(dbz5, dbz5<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')
ax.contourf(xh6, yh, np.ma.masked_array(dbz6, dbz6<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')
ax.contourf(xh7, yh, np.ma.masked_array(dbz7, dbz7<20), levels=dbz_levs, vmin=0, vmax=70, cmap='HomeyerRainbow')

cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_ticks(np.linspace(0,70,8))
cb.set_label('Radar reflectivity (dBZ)', fontsize=10)

ax.contour(xh1, yh, hail1, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh2, yh, hail2, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh3, yh, hail3, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh4, yh, hail4, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh5, yh, hail5, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh6, yh, hail6, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh7, yh, hail7, levels=levs, colors=cols, linewidths=lws)
ax.set_xlim([0,450])
ax.set_ylim([60,180])
ax.set_xlabel('Translated x (km)', fontsize=10)
ax.set_ylabel('y (km)', fontsize=10)
ax.set_title('Brooks test run #1: Translated accumulated solid precip swaths (0-7 h)', fontsize=12)
l1, = ax.plot([-2,-1], [-2,-1], 'gray', linewidth=1)
l2, = ax.plot([-2,-1], [-2,-1], 'k', linewidth=1)
ax.legend(handles=[l1,l2], labels=['0.1 mm','1 mm'], loc='lower right', fontsize=9)

# plt.savefig(fp+'dbz_hail_swath_translated.png', dpi=300)



plt.show()




# I think i have to interpolate dear god


#%% Check tendency nudging

S0 = -0.015
t1 = 1800
t2 = 2400
t3 = 2401
t4 = 3000
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
ax.set_ylim([-15,15])

plt.show()











