# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 15:05:24 2026

@author: mschne28
"""

from CM1utils import *

#%% Plot simulations


fp = 'D:/brooks/rerun2/'
# figstr = 'ERA5-1_final'

ds = xr.open_dataset(fp+'cm1out_000029.nc')
time = ds['time'][:].values[0]
xh = ds['xh'][:].values
yh = ds['yh'][:].values
zh = ds['zh'][:].values
iz3 = np.argmin(abs(zh-3))
iz1 = np.argmin(abs(zh-1))

umove = ds['umove'].values[0]
vmove = ds['vmove'].values[0]
dbz = ds.sel(time=time, zh=zh[0])['dbz'][:].values
thpert = ds.sel(time=time, zh=zh[0])['th'][:].values - ds.sel(time=time, zh=zh[0])['th0'][:].values
winterp = ds.sel(time=time, zh=slice(zh[0],zh[iz3+1]))['winterp'][:].values
zvort = ds.sel(time=time, zh=slice(zh[0],zh[iz1+1]))['zvort'][:].values
# uinterp = ds.sel(time=time, zh=zh[0])['uinterp'][:].values + umove
# vinterp = ds.sel(time=time, zh=zh[0])['vinterp'][:].values + vmove

ds.close()



iz = np.where(zh>=3)[0][0]
wmax = np.max(winterp[0:iz,:,:], axis=0)

iz = np.where(zh>=1)[0][0]
vortmax = np.max(zvort[0:iz,:,:], axis=0)

# xh = xh + 2*time/1000


time = time / np.timedelta64(1, 's')


figsave = False



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
if figsave:
    plt.savefig(fp+f"figs/dbz_{time/60:.0f}min.png", dpi=300)


# plt.show()



# figsave = False



# fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
# plot_contourf(xh, yh, wmax, 'w', ax, levels=np.linspace(-15,15,31), datalims=[-15,15],
#               cmap='balance', cbticks=np.linspace(-15,15,11), cbfs=12)
# ax.contour(xh, yh, vortmax, levels=[0.025,0.05], colors=['b','dodgerblue'], linewidths=[1,1.5])
# ax.set_xlim([-100,100])
# ax.set_ylim([-100,100])
# ax.set_title(f"Max 0-3 km w, max 0-1 km \u03B6 (t = {time/60:.0f} min)", fontsize=12)
# ax.set_xlabel('x (km)', fontsize=12)
# ax.set_ylabel('y (km)', fontsize=12)

# l1,= ax.plot([-150,-149], [-150,-149], color='b', linestyle='-', linewidth=1)
# l2,= ax.plot([-150,-149], [-150,-149], color='dodgerblue', linestyle='-', linewidth=1.5)
# # l3,= ax.plot([-150,-149], [-150,-149], color='k', linestyle='-', linewidth=1)
# ax.legend(handles=[l1,l2], labels=["\u03B6=0.025 s$^{-1}$", "\u03B6=0.05 s$^{-1}$"],
#           fontsize=10, loc='lower right')
# if figsave:
#     plt.savefig(fp+f"figs/wmax_{time/60:.0f}min.png", dpi=300)



# qix = 40
# fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
# plot_contourf(xh, yh, thpert, 'thpert', ax, levels=np.linspace(-10,10,21), datalims=[-10,10],
#               cmap='balance', cbticks=np.linspace(-10,10,11), cbfs=12)
# # ax.quiver(xh[::qix], yh[::qix], uinterp[::qix,::qix], vinterp[::qix,::qix], color='k', scale=200, width=0.003, pivot='middle')
# ax.set_xlim([-100,100])
# ax.set_ylim([-100,100])
# ax.set_title(f"Surface \u03B8', surface winds (t = {time/60:.0f} min)", fontsize=12)
# ax.set_xlabel('x (km)', fontsize=12)
# ax.set_ylabel('y (km)', fontsize=12)
# # if figsave:
# #     plt.savefig(fp+f"figs/thpert_{time/60:.0f}min.png", dpi=300)



# plt.show()


i_wmax = np.where(wmax == np.max(wmax))
i_vmax = np.where(vortmax == np.max(vortmax))

x_mc = np.mean([xh[i_wmax[1][0]], xh[i_vmax[1][0]]])
y_mc = np.mean([yh[i_wmax[0][0]], yh[i_vmax[0][0]]])

u_trans = x_mc*1000/time
v_trans = y_mc*1000/time

umove_new = u_trans + umove
vmove_new = v_trans + vmove


print(f"Time = {time/60:.0f} min")
# print(f"Old umove = {umove:.1f} m/s || Old vmove = {vmove:.1f} m/s")
print(f"New umove = {umove_new:.1f} m/s")
print(f"New vmove = {vmove_new:.1f} m/s")
print("-----")








#%% Check stats time series

fp = 'D:/brooks/rerun1/'

ds = nc.Dataset(fp+'cm1out_stats.nc')
time = ds.variables['time'][:].data + 120
vortsfc = ds.variables['vortsfc'][:].data
vort1km = ds.variables['vort1km'][:].data
vort3km = ds.variables['vort3km'][:].data
vort5km = ds.variables['vort5km'][:].data
wmax1000 = ds.variables['wmax1000'][:].data
wmax2500 = ds.variables['wmax2500'][:].data
wmax5000 = ds.variables['wmax5000'][:].data
swspmax = ds.variables['swspmax'][:].data
wsp1max = ds.variables['wsp1max'][:].data
wsp80max = ds.variables['wsp0max'][:].data
wsp1min = ds.variables['wsp1min'][:].data
wsp0min = ds.variables['wsp0min'][:].data
dmimax = np.max(np.asarray([ds.variables['smaxdmi1'][:].data, ds.variables['smaxdmi2'][:].data, ds.variables['smaxdmi3'][:].data, ds.variables['smaxdmi4'][:].data]), axis=0)
dhlmax = np.max(np.asarray([ds.variables['smaxdhl1'][:].data, ds.variables['smaxdhl2'][:].data, ds.variables['smaxdhl3'][:].data, ds.variables['smaxdhl4'][:].data]), axis=0)
# dmimax = np.max(np.asarray([ds.variables['maxdmi1'][:].data, ds.variables['maxdmi2'][:].data, ds.variables['maxdmi3'][:].data, ds.variables['maxdmi4'][:].data]), axis=0)
# dhlmax = np.max(np.asarray([ds.variables['maxdhl1'][:].data, ds.variables['maxdhl2'][:].data, ds.variables['maxdhl3'][:].data, ds.variables['maxdhl4'][:].data]), axis=0)
zwmax = ds.variables['zwmax'][:].data
zvortmax = ds.variables['zvortmax'][:].data
ds.close()


figsave = False




fig,ax = plt.subplots(1, 1, figsize=(10,4), layout='constrained')
l1, = ax.plot(time/60, movmean(vort5km,5), 'k', linewidth=2)
l2, = ax.plot(time/60, movmean(vort3km,5), 'b', linewidth=2)
l3, = ax.plot(time/60, movmean(vort1km,5), 'deepskyblue', linewidth=2)
l4, = ax.plot(time/60, movmean(vortsfc,5), 'r', linewidth=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax.xaxis.set_major_locator(MultipleLocator(30))
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(0.04))
ax.yaxis.set_minor_locator(MultipleLocator(0.01))
ax.set_xlim([0,420])
ax.set_ylim([0,0.2])
ax.legend(handles=[l4,l3,l2,l1], labels=["\u03B6(sfc)", "\u03B6(1 km)", "\u03B6(3 km)", "\u03B6(5 km)"], loc='lower right', fontsize=12, ncols=2)
ax.set_title("Domain-maximum \u03B6", fontsize=14)
ax.set_xlabel('Time (min)', fontsize=12)
ax.set_ylabel("\u03B6 (1/s)", fontsize=12)
if figsave:
    plt.savefig(fp+f"figs/vortmax_timeseries.png", dpi=300)


fig,ax = plt.subplots(1, 1, figsize=(10,4), layout='constrained')
l1, = ax.plot(time/60, movmean(wmax1000,5), 'deepskyblue', linewidth=2)
l2, = ax.plot(time/60, movmean(wmax2500,5), 'b', linewidth=2)
l3, = ax.plot(time/60, movmean(wmax5000,5), 'k', linewidth=2)
ax.grid(visible=True, which='major', color='dimgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax.xaxis.set_major_locator(MultipleLocator(30))
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.set_xlim([0,420])
ax.set_ylim([0,30])
ax.legend(handles=[l1,l2,l3], labels=['w(1 km)', 'w(2.5 km)', 'w(5 km)'], loc='lower right', fontsize=12)
ax.set_title("Domain-maximum w", fontsize=14)
ax.set_xlabel('Time (min)', fontsize=12)
ax.set_ylabel('w (m/s)', fontsize=12)
if figsave:
    plt.savefig(fp+f"figs/wmax_timeseries.png", dpi=300)


fig,ax = plt.subplots(1, 1, figsize=(10,4), layout='constrained')
l1, = ax.plot(time/60, movmean(swspmax,5), 'deepskyblue', linewidth=2)
l2, = ax.plot(time/60, movmean(wsp80max,5), 'b', linewidth=2)
l3, = ax.plot(time/60, movmean(wsp1max,5), 'k', linewidth=2)
ax.plot([0,1000], [25.7,25.7], '--k', linewidth=1)
ax.grid(visible=True, which='major', color='dimgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax.xaxis.set_major_locator(MultipleLocator(30))
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.set_xlim([0,420])
ax.set_ylim([0,40])
ax.legend(handles=[l1,l2,l3], labels=['wsp(10 m)', 'wsp(80 m)', 'wsp(1 km)'], loc='lower right', fontsize=12)
ax.set_title("Domain-maximum wspd", fontsize=14)
ax.set_xlabel('Time (min)', fontsize=12)
ax.set_ylabel('wspd (m/s)', fontsize=12)


fig,ax = plt.subplots(1, 1, figsize=(10,4), layout='constrained')
l1, = ax.plot(time/60, movmean(zwmax/1000,10), 'dimgray', linewidth=2)
l2, = ax.plot(time/60, movmean(zvortmax/1000,10), 'k', linewidth=2)
axx = ax.twinx()
l3, = axx.plot(time/60, movmean(dmimax*100,5), 'r', linewidth=2)
l4, = axx.plot(time/60, movmean(dhlmax*100,5), 'b', linewidth=2)
ax.grid(visible=True, which='major', color='dimgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax.xaxis.set_major_locator(MultipleLocator(30))
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.set_xlim([0,420])
ax.set_ylim([0,10])
axx.set_ylim([0,10])
ax.legend(handles=[l1,l2,l3,l4], labels=["z(max w)", "z(max \u03B6)", "D$_{ice}$", "D$_{hail}$"], loc='upper right', fontsize=12, ncols=2)
ax.set_title("Height of max w/\u03B6, max sfc ice and hail diameter", fontsize=14)
ax.set_xlabel('Time (min)', fontsize=12)
ax.set_ylabel('Height (km)', fontsize=12)
axx.set_ylabel('Diameter (cm)', fontsize=12)

plt.show()



#%% Compare stat time series in test simulations (just for testing)

fp = "D:/brooks/"

ds = nc.Dataset(fp+"era5-1_125m_test2/cm1out_stats.nc")
time = ds.variables['time'][:].data + 120
vortsfc_1a = ds.variables['vortsfc'][:].data
vort1km_1a = ds.variables['vort1km'][:].data
vort3km_1a = ds.variables['vort3km'][:].data
vort5km_1a = ds.variables['vort5km'][:].data
wmax1000_1a = ds.variables['wmax1000'][:].data
wmax2500_1a = ds.variables['wmax2500'][:].data
wmax5000_1a = ds.variables['wmax5000'][:].data
ds.close()


ds = nc.Dataset(fp+"era5-1_125m_test3/cm1out_stats.nc")
vortsfc_2a = ds.variables['vortsfc'][:].data
vort1km_2a = ds.variables['vort1km'][:].data
vort3km_2a = ds.variables['vort3km'][:].data
vort5km_2a = ds.variables['vort5km'][:].data
wmax1000_2a = ds.variables['wmax1000'][:].data
wmax2500_2a = ds.variables['wmax2500'][:].data
wmax5000_2a = ds.variables['wmax5000'][:].data
ds.close()


ds = nc.Dataset(fp+"era5-1_125m_test4/cm1out_stats.nc")
vortsfc_3a = ds.variables['vortsfc'][:].data
vort1km_3a = ds.variables['vort1km'][:].data
vort3km_3a = ds.variables['vort3km'][:].data
vort5km_3a = ds.variables['vort5km'][:].data
wmax1000_3a = ds.variables['wmax1000'][:].data
wmax2500_3a = ds.variables['wmax2500'][:].data
wmax5000_3a = ds.variables['wmax5000'][:].data
ds.close()


ds = nc.Dataset(fp+"era5-1_125m_test5/cm1out_stats.nc")
vortsfc_4a = ds.variables['vortsfc'][:].data
vort1km_4a = ds.variables['vort1km'][:].data
vort3km_4a = ds.variables['vort3km'][:].data
vort5km_4a = ds.variables['vort5km'][:].data
wmax1000_4a = ds.variables['wmax1000'][:].data
wmax2500_4a = ds.variables['wmax2500'][:].data
wmax5000_4a = ds.variables['wmax5000'][:].data
ds.close()


ds = nc.Dataset(fp+"era5-1_125m_test6/cm1out_stats.nc")
vortsfc_5a = ds.variables['vortsfc'][:].data
vort1km_5a = ds.variables['vort1km'][:].data
vort3km_5a = ds.variables['vort3km'][:].data
vort5km_5a = ds.variables['vort5km'][:].data
wmax1000_5a = ds.variables['wmax1000'][:].data
wmax2500_5a = ds.variables['wmax2500'][:].data
wmax5000_5a = ds.variables['wmax5000'][:].data
ds.close()


# ds = nc.Dataset(fp+"era5-1_125m_test7/cm1out_stats.nc")
# vortsfc_6a = ds.variables['vortsfc'][:].data
# vort1km_6a = ds.variables['vort1km'][:].data
# vort3km_6a = ds.variables['vort3km'][:].data
# vort5km_6a = ds.variables['vort5km'][:].data
# wmax1000_6a = ds.variables['wmax1000'][:].data
# wmax2500_6a = ds.variables['wmax2500'][:].data
# wmax5000_6a = ds.variables['wmax5000'][:].data
# ds.close()



# ds = nc.Dataset(fp+"era5-2_125m_test1/cm1out_stats.nc")
# time = ds.variables['time'][:].data + 120
# vortsfc_1b = ds.variables['vortsfc'][:].data
# vort1km_1b = ds.variables['vort1km'][:].data
# vort3km_1b = ds.variables['vort3km'][:].data
# vort5km_1b = ds.variables['vort5km'][:].data
# wmax1000_1b = ds.variables['wmax1000'][:].data
# wmax2500_1b = ds.variables['wmax2500'][:].data
# wmax5000_1b = ds.variables['wmax5000'][:].data
# ds.close()



# ds = nc.Dataset(fp+"era5-3_125m_test1/cm1out_stats.nc")
# time = ds.variables['time'][:].data + 120
# vortsfc_1c = ds.variables['vortsfc'][:].data
# vort1km_1c = ds.variables['vort1km'][:].data
# vort3km_1c = ds.variables['vort3km'][:].data
# vort5km_1c = ds.variables['vort5km'][:].data
# wmax1000_1c = ds.variables['wmax1000'][:].data
# wmax2500_1c = ds.variables['wmax2500'][:].data
# wmax5000_1c = ds.variables['wmax5000'][:].data
# ds.close()



figsave = False


### Compare each variable across simulations, separated by variable



fig,ax = plt.subplots(4, 1, figsize=(9,10), sharex=True, layout='constrained')

# l1, = ax[0].plot(time/60, movmean(vort5km_1a,5), 'k', linewidth=2)
# l2, = ax[0].plot(time/60, movmean(vort5km_2a,5), 'b', linewidth=2)
l3, = ax[0].plot(time/60, movmean(vort5km_3a,5), 'deepskyblue', linewidth=2)
l4, = ax[0].plot(time/60, movmean(vort5km_4a,5), 'darkviolet', linewidth=2)
l5, = ax[0].plot(time/60, movmean(vort5km_5a,5), 'violet', linewidth=2)
# l6, = ax[0].plot(time/60, movmean(vort5km_6a,5), 'magenta', linewidth=2)
# hh = [l1,l2,l3,l4,l5]
# ll = ["\u03B8'=9 K", "\u03B8'=4.5 K", "\u03B8'=3 K", "\u03B8'=3 K (25km)", "\u03B8'=3 K (30km)"]
# hh = [l1,l2,l3]
# ll = ["\u03B8'=9 K", "\u03B8'=4.5 K", "\u03B8'=3 K"]
hh = [l3,l4,l5]
ll = ["\u03B8'=3 K", "\u03B8'=3 K (25km)", "\u03B8'=3 K (30km)"]
ax[0].legend(handles=hh, labels=ll, loc='upper right', fontsize=11, ncols=2)
ax[0].set_title("Domain-maximum \u03B6(5 km)", fontsize=14)
ax[0].set_ylabel("\u03B6 (1/s)", fontsize=12)
ax[0].set_ylim([0,0.2])

# ax[1].plot(time/60, movmean(vort3km_1a,5), 'k', linewidth=2)
# ax[1].plot(time/60, movmean(vort3km_2a,5), 'b', linewidth=2)
ax[1].plot(time/60, movmean(vort3km_3a,5), 'deepskyblue', linewidth=2)
ax[1].plot(time/60, movmean(vort3km_4a,5), 'darkviolet', linewidth=2)
ax[1].plot(time/60, movmean(vort3km_5a,5), 'violet', linewidth=2)
# ax[1].plot(time/60, movmean(vort3km_6a,5), 'magenta', linewidth=2)
ax[1].set_title("Domain-maximum \u03B6 (3 km)", fontsize=14)
ax[1].set_ylabel("\u03B6 (1/s)", fontsize=12)
ax[1].set_ylim([0,0.2])

# ax[2].plot(time/60, movmean(vort1km_1a,5), 'k', linewidth=2)
# ax[2].plot(time/60, movmean(vort1km_2a,5), 'b', linewidth=2)
ax[2].plot(time/60, movmean(vort1km_3a,5), 'deepskyblue', linewidth=2)
ax[2].plot(time/60, movmean(vort1km_4a,5), 'darkviolet', linewidth=2)
ax[2].plot(time/60, movmean(vort1km_5a,5), 'violet', linewidth=2)
# ax[2].plot(time/60, movmean(vort1km_6a,5), 'magenta', linewidth=2)
ax[2].set_title("Domain-maximum \u03B6(1 km)", fontsize=14)
ax[2].set_ylabel("\u03B6 (1/s)", fontsize=12)
ax[2].set_ylim([0,0.2])

# ax[3].plot(time/60, movmean(vortsfc_1a,5), 'k', linewidth=2)
# ax[3].plot(time/60, movmean(vortsfc_2a,5), 'b', linewidth=2)
ax[3].plot(time/60, movmean(vortsfc_3a,5), 'deepskyblue', linewidth=2)
ax[3].plot(time/60, movmean(vortsfc_4a,5), 'darkviolet', linewidth=2)
ax[3].plot(time/60, movmean(vortsfc_5a,5), 'violet', linewidth=2)
# ax[3].plot(time/60, movmean(vortsfc_6a,5), 'magenta', linewidth=2)
ax[3].set_title("Domain-maximum \u03B6(sfc)", fontsize=14)
ax[3].set_xlabel('Time (min)', fontsize=12)
ax[3].set_ylabel("\u03B6 (1/s)", fontsize=12)
ax[3].set_xlim([0,480])
ax[3].set_ylim([0,0.28])

for n in range(len(ax)):
    ax[n].grid(visible=True, which='major', color='gray', linestyle='-')
    ax[n].grid(visible=True, which='minor', color='lightgray', linestyle='-')
    ax[n].xaxis.set_major_locator(MultipleLocator(30))
    ax[n].xaxis.set_minor_locator(MultipleLocator(10))
    ax[n].yaxis.set_major_locator(MultipleLocator(0.04))
    ax[n].yaxis.set_minor_locator(MultipleLocator(0.02))
    
if figsave:
    plt.savefig(fp+f"figs/vortmax_timeseries_ERA5-1_compare_1.png", dpi=300)






fig,ax = plt.subplots(3, 1, figsize=(9,9), sharex=True, layout='constrained')

# l1, = ax[0].plot(time/60, movmean(wmax5000_1a,5), 'k', linewidth=2)
# l2, = ax[0].plot(time/60, movmean(wmax5000_2a,5), 'b', linewidth=2)
l3, = ax[0].plot(time/60, movmean(wmax5000_3a,5), 'deepskyblue', linewidth=2)
l4, = ax[0].plot(time/60, movmean(wmax5000_4a,5), 'darkviolet', linewidth=2)
l5, = ax[0].plot(time/60, movmean(wmax5000_5a,5), 'violet', linewidth=2)
ax[0].legend(handles=hh, labels=ll, loc='lower right', fontsize=11, ncols=2)
ax[0].set_title("Domain-maximum w(5 km)", fontsize=14)
ax[0].set_ylabel("w (m/s)", fontsize=12)
ax[0].set_ylim([0,25])

# ax[1].plot(time/60, movmean(wmax2500_1a,5), 'k', linewidth=2)
# ax[1].plot(time/60, movmean(wmax2500_2a,5), 'b', linewidth=2)
ax[1].plot(time/60, movmean(wmax2500_3a,5), 'deepskyblue', linewidth=2)
ax[1].plot(time/60, movmean(wmax2500_4a,5), 'darkviolet', linewidth=2)
ax[1].plot(time/60, movmean(wmax2500_5a,5), 'violet', linewidth=2)
ax[1].set_title("Domain-maximum w(2.5 km)", fontsize=14)
ax[1].set_ylabel("w (m/s)", fontsize=12)
ax[1].set_ylim([0,25])

# ax[2].plot(time/60, movmean(wmax1000_1a,5), 'k', linewidth=2)
# ax[2].plot(time/60, movmean(wmax1000_2a,5), 'b', linewidth=2)
ax[2].plot(time/60, movmean(wmax1000_3a,5), 'deepskyblue', linewidth=2)
ax[2].plot(time/60, movmean(wmax1000_4a,5), 'darkviolet', linewidth=2)
ax[2].plot(time/60, movmean(wmax1000_5a,5), 'violet', linewidth=2)
ax[2].set_title("Domain-maximum w(1 km)", fontsize=14)
ax[2].set_xlabel('Time (min)', fontsize=12)
ax[2].set_ylabel("w (m/s)", fontsize=12)
ax[2].set_xlim([0,480])
ax[2].set_ylim([0,25])

for n in range(len(ax)):
    ax[n].grid(visible=True, which='major', color='gray', linestyle='-')
    ax[n].grid(visible=True, which='minor', color='lightgray', linestyle='-')
    ax[n].xaxis.set_major_locator(MultipleLocator(30))
    ax[n].xaxis.set_minor_locator(MultipleLocator(10))
    ax[n].yaxis.set_major_locator(MultipleLocator(5))
    ax[n].yaxis.set_minor_locator(MultipleLocator(1))
    
if figsave:
    plt.savefig(fp+f"figs/wmax_timeseries_ERA5-1_compare_1.png", dpi=300)





### Compare variables within each simulation, separated by simulation

fig,ax = plt.subplots(3, 1, figsize=(9,9), sharex=True, sharey=True, layout='constrained')

l1, = ax[0].plot(time/60, movmean(vort5km_1a,5), 'k', linewidth=2)
l2, = ax[0].plot(time/60, movmean(vort3km_1a,5), 'b', linewidth=2)
l3, = ax[0].plot(time/60, movmean(vort1km_1a,5), 'deepskyblue', linewidth=2)
l4, = ax[0].plot(time/60, movmean(vortsfc_1a,5), 'magenta', linewidth=1.5)
ax[0].set_title("Domain-maximum \u03B6 (\u03B8'=9 K)", fontsize=14)
ax[0].set_ylabel("\u03B6 (1/s)", fontsize=12)
ax[0].set_xlim([0,480])
ax[0].set_ylim([0,0.28])
ax[0].legend(handles=[l4,l3,l2,l1], labels=["\u03B6(sfc)", "\u03B6(1 km)", "\u03B6(3 km)", "\u03B6(5 km)"], loc='upper right', fontsize=11, ncols=2)

ax[1].plot(time/60, movmean(vort5km_2a,5), 'k', linewidth=2)
ax[1].plot(time/60, movmean(vort3km_2a,5), 'b', linewidth=2)
ax[1].plot(time/60, movmean(vort1km_2a,5), 'deepskyblue', linewidth=2)
ax[1].plot(time/60, movmean(vortsfc_2a,5), 'magenta', linewidth=1.5)
ax[1].set_title("Domain-maximum \u03B6 (\u03B8'=4.5 K)", fontsize=14)
ax[1].set_ylabel("\u03B6 (1/s)", fontsize=12)

ax[2].plot(time/60, movmean(vort5km_3a,5), 'k', linewidth=2)
ax[2].plot(time/60, movmean(vort3km_3a,5), 'b', linewidth=2)
ax[2].plot(time/60, movmean(vort1km_3a,5), 'deepskyblue', linewidth=2)
ax[2].plot(time/60, movmean(vortsfc_3a,5), 'magenta', linewidth=1.5)
ax[2].set_title("Domain-maximum \u03B6 (\u03B8'=3 K)", fontsize=14)
ax[2].set_xlabel('Time (min)', fontsize=12)
ax[2].set_ylabel("\u03B6 (1/s)", fontsize=12)

for n in range(len(ax)):
    ax[n].grid(visible=True, which='major', color='gray', linestyle='-')
    ax[n].grid(visible=True, which='minor', color='lightgray', linestyle='-')
    ax[n].xaxis.set_major_locator(MultipleLocator(30))
    ax[n].xaxis.set_minor_locator(MultipleLocator(10))
    ax[n].yaxis.set_major_locator(MultipleLocator(0.04))
    ax[n].yaxis.set_minor_locator(MultipleLocator(0.02))
    
if figsave:
    plt.savefig(fp+f"figs/vortmax_timeseries_ERA5-1_compare_2.png", dpi=300)




fig,ax = plt.subplots(3, 1, figsize=(9,9), sharex=True, sharey=True, layout='constrained')

l1, = ax[0].plot(time/60, movmean(wmax5000_1a,5), 'k', linewidth=2)
l2, = ax[0].plot(time/60, movmean(wmax2500_1a,5), 'b', linewidth=2)
l3, = ax[0].plot(time/60, movmean(wmax1000_1a,5), 'deepskyblue', linewidth=2)
ax[0].set_title("Domain-maximum w (\u03B8'=9 K)", fontsize=14)
ax[0].set_ylabel("w (m/s)", fontsize=12)
ax[0].set_xlim([0,480])
ax[0].set_ylim([0,25])
ax[0].legend(handles=[l3,l2,l1], labels=["w(1 km)", "w(2.5 km)", "w(3 km)"], loc='lower right', fontsize=11)

ax[1].plot(time/60, movmean(wmax5000_2a,5), 'k', linewidth=2)
ax[1].plot(time/60, movmean(wmax2500_2a,5), 'b', linewidth=2)
ax[1].plot(time/60, movmean(wmax1000_2a,5), 'deepskyblue', linewidth=2)
ax[1].set_title("Domain-maximum w (\u03B8'=4.5 K)", fontsize=14)
ax[1].set_ylabel("w (m/s)", fontsize=12)

ax[2].plot(time/60, movmean(wmax5000_3a,5), 'k', linewidth=2)
ax[2].plot(time/60, movmean(wmax2500_3a,5), 'b', linewidth=2)
ax[2].plot(time/60, movmean(wmax1000_3a,5), 'deepskyblue', linewidth=2)
ax[2].set_title("Domain-maximum w (\u03B8'=3 K)", fontsize=14)
ax[2].set_xlabel('Time (min)', fontsize=12)
ax[2].set_ylabel("w (m/s)", fontsize=12)

for n in range(len(ax)):
    ax[n].grid(visible=True, which='major', color='gray', linestyle='-')
    ax[n].grid(visible=True, which='minor', color='lightgray', linestyle='-')
    ax[n].xaxis.set_major_locator(MultipleLocator(30))
    ax[n].xaxis.set_minor_locator(MultipleLocator(10))
    ax[n].yaxis.set_major_locator(MultipleLocator(5))
    ax[n].yaxis.set_minor_locator(MultipleLocator(1))
    
if figsave:
    plt.savefig(fp+f"figs/wmax_timeseries_ERA5-1_compare_2.png", dpi=300)





#%% Compare stat time series with bubble spacing

fig,ax = plt.subplots(3, 1, figsize=(9,9), sharex=True, sharey=True, layout='constrained')

l1, = ax[0].plot(time/60, movmean(vort5km_3a,5), 'k', linewidth=2)
l2, = ax[0].plot(time/60, movmean(vort3km_3a,5), 'b', linewidth=2)
l3, = ax[0].plot(time/60, movmean(vort1km_3a,5), 'deepskyblue', linewidth=2)
l4, = ax[0].plot(time/60, movmean(vortsfc_3a,5), 'magenta', linewidth=1.5)
ax[0].set_title("Domain-maximum \u03B6 (\u03B8'=3 K, dy=20 km)", fontsize=14)
ax[0].set_ylabel("\u03B6 (1/s)", fontsize=12)
ax[0].set_xlim([0,480])
ax[0].set_ylim([0,0.2])
ax[0].legend(handles=[l4,l3,l2,l1], labels=["\u03B6(sfc)", "\u03B6(1 km)", "\u03B6(3 km)", "\u03B6(5 km)"], loc='lower right', fontsize=11, ncols=2)

ax[1].plot(time/60, movmean(vort5km_4a,5), 'k', linewidth=2)
ax[1].plot(time/60, movmean(vort3km_4a,5), 'b', linewidth=2)
ax[1].plot(time/60, movmean(vort1km_4a,5), 'deepskyblue', linewidth=2)
ax[1].plot(time/60, movmean(vortsfc_4a,5), 'magenta', linewidth=1.5)
ax[1].set_title("Domain-maximum \u03B6 (\u03B8'=3 K, dy=25 km)", fontsize=14)
ax[1].set_ylabel("\u03B6 (1/s)", fontsize=12)

ax[2].plot(time/60, movmean(vort5km_5a,5), 'k', linewidth=2)
ax[2].plot(time/60, movmean(vort3km_5a,5), 'b', linewidth=2)
ax[2].plot(time/60, movmean(vort1km_5a,5), 'deepskyblue', linewidth=2)
ax[2].plot(time/60, movmean(vortsfc_5a,5), 'magenta', linewidth=1.5)
ax[2].set_title("Domain-maximum \u03B6 (\u03B8'=3 K, dy=30 km)", fontsize=14)
ax[2].set_xlabel('Time (min)', fontsize=12)
ax[2].set_ylabel("\u03B6 (1/s)", fontsize=12)

for n in range(len(ax)):
    ax[n].grid(visible=True, which='major', color='gray', linestyle='-')
    ax[n].grid(visible=True, which='minor', color='lightgray', linestyle='-')
    ax[n].xaxis.set_major_locator(MultipleLocator(30))
    ax[n].xaxis.set_minor_locator(MultipleLocator(10))
    ax[n].yaxis.set_major_locator(MultipleLocator(0.04))
    ax[n].yaxis.set_minor_locator(MultipleLocator(0.02))



fig,ax = plt.subplots(3, 1, figsize=(9,9), sharex=True, sharey=True, layout='constrained')

l1, = ax[0].plot(time/60, movmean(wmax5000_3a,5), 'k', linewidth=2)
l2, = ax[0].plot(time/60, movmean(wmax2500_3a,5), 'b', linewidth=2)
l3, = ax[0].plot(time/60, movmean(wmax1000_3a,5), 'deepskyblue', linewidth=2)
ax[0].set_title("Domain-maximum w (\u03B8'=3 K, dy=20 km)", fontsize=14)
ax[0].set_ylabel("w (m/s)", fontsize=12)
ax[0].set_xlim([0,480])
ax[0].set_ylim([0,25])
ax[0].legend(handles=[l3,l2,l1], labels=["w(1 km)", "w(2.5 km)", "w(5 km)"], loc='lower right', fontsize=11)

ax[1].plot(time/60, movmean(wmax5000_4a,5), 'k', linewidth=2)
ax[1].plot(time/60, movmean(wmax2500_4a,5), 'b', linewidth=2)
ax[1].plot(time/60, movmean(wmax1000_4a,5), 'deepskyblue', linewidth=2)
ax[1].set_title("Domain-maximum w (\u03B8'=3 K, dy=25 km)", fontsize=14)
ax[1].set_ylabel("w (m/s)", fontsize=12)

ax[2].plot(time/60, movmean(wmax5000_5a,5), 'k', linewidth=2)
ax[2].plot(time/60, movmean(wmax2500_5a,5), 'b', linewidth=2)
ax[2].plot(time/60, movmean(wmax1000_5a,5), 'deepskyblue', linewidth=2)
ax[2].set_title("Domain-maximum w (\u03B8'=3 K, dy=30 km)", fontsize=14)
ax[2].set_xlabel('Time (min)', fontsize=12)
ax[2].set_ylabel("w (m/s)", fontsize=12)

for n in range(len(ax)):
    ax[n].grid(visible=True, which='major', color='gray', linestyle='-')
    ax[n].grid(visible=True, which='minor', color='lightgray', linestyle='-')
    ax[n].xaxis.set_major_locator(MultipleLocator(30))
    ax[n].xaxis.set_minor_locator(MultipleLocator(10))
    ax[n].yaxis.set_major_locator(MultipleLocator(5))
    ax[n].yaxis.set_minor_locator(MultipleLocator(1))


plt.show()



#%% Compare between soundings

# fig,ax = plt.subplots(4, 1, figsize=(9,10), sharex=True, layout='constrained')

# l1, = ax[0].plot(time/60, movmean(vort5km_5a,5), 'k', linewidth=2)
# l2, = ax[0].plot(time/60, movmean(vort5km_1b,5), 'b', linewidth=2)
# l3, = ax[0].plot(time/60, movmean(vort5km_1c,5), 'r', linewidth=2)
# hh = [l1,l2,l3]
# ll = ["ERA5-1", "ERA5-2", "ERA5-3"]
# ax[0].legend(handles=hh, labels=ll, loc='upper right', fontsize=11, ncols=2)
# ax[0].set_title("Domain-maximum \u03B6(5 km)", fontsize=14)
# ax[0].set_ylabel("\u03B6 (1/s)", fontsize=12)
# ax[0].set_ylim([0,0.2])

# ax[1].plot(time/60, movmean(vort3km_5a,5), 'k', linewidth=2)
# ax[1].plot(time/60, movmean(vort3km_1b,5), 'b', linewidth=2)
# ax[1].plot(time/60, movmean(vort3km_1c,5), 'r', linewidth=2)
# ax[1].set_title("Domain-maximum \u03B6 (3 km)", fontsize=14)
# ax[1].set_ylabel("\u03B6 (1/s)", fontsize=12)
# ax[1].set_ylim([0,0.2])

# ax[2].plot(time/60, movmean(vort1km_5a,5), 'k', linewidth=2)
# ax[2].plot(time/60, movmean(vort1km_1b,5), 'b', linewidth=2)
# ax[2].plot(time/60, movmean(vort1km_1c,5), 'r', linewidth=2)
# ax[2].set_title("Domain-maximum \u03B6(1 km)", fontsize=14)
# ax[2].set_ylabel("\u03B6 (1/s)", fontsize=12)
# ax[2].set_ylim([0,0.2])

# ax[3].plot(time/60, movmean(vortsfc_5a,5), 'k', linewidth=2)
# ax[3].plot(time/60, movmean(vortsfc_1b,5), 'b', linewidth=2)
# ax[3].plot(time/60, movmean(vortsfc_1c,5), 'r', linewidth=2)
# ax[3].set_title("Domain-maximum \u03B6(sfc)", fontsize=14)
# ax[3].set_xlabel('Time (min)', fontsize=12)
# ax[3].set_ylabel("\u03B6 (1/s)", fontsize=12)
# ax[3].set_xlim([0,480])
# ax[3].set_ylim([0,0.28])

# for n in range(len(ax)):
#     ax[n].grid(visible=True, which='major', color='gray', linestyle='-')
#     ax[n].grid(visible=True, which='minor', color='lightgray', linestyle='-')
#     ax[n].xaxis.set_major_locator(MultipleLocator(30))
#     ax[n].xaxis.set_minor_locator(MultipleLocator(10))
#     ax[n].yaxis.set_major_locator(MultipleLocator(0.04))
#     ax[n].yaxis.set_minor_locator(MultipleLocator(0.02))
    
# if figsave:
#     plt.savefig(fp+'vortmax_timeseries_ERA5_compare_3.png', dpi=300)



# fig,ax = plt.subplots(3, 1, figsize=(9,9), sharex=True, layout='constrained')

# l1, = ax[0].plot(time/60, movmean(wmax5000_5a,5), 'k', linewidth=2)
# l2, = ax[0].plot(time/60, movmean(wmax5000_1b,5), 'b', linewidth=2)
# l3, = ax[0].plot(time/60, movmean(wmax5000_1c,5), 'r', linewidth=2)
# hh = [l1,l2,l3]
# ll = ["ERA5-1", "ERA5-2", "ERA5-3"]
# ax[0].legend(handles=hh, labels=ll, loc='lower right', fontsize=11, ncols=2)
# ax[0].set_title("Domain-maximum w(5 km)", fontsize=14)
# ax[0].set_ylabel("w (m/s)", fontsize=12)
# ax[0].set_ylim([0,25])

# ax[1].plot(time/60, movmean(wmax2500_5a,5), 'k', linewidth=2)
# ax[1].plot(time/60, movmean(wmax2500_1b,5), 'b', linewidth=2)
# ax[1].plot(time/60, movmean(wmax2500_1c,5), 'r', linewidth=2)
# ax[1].set_title("Domain-maximum w(2.5 km)", fontsize=14)
# ax[1].set_ylabel("w (m/s)", fontsize=12)
# ax[1].set_ylim([0,25])

# ax[2].plot(time/60, movmean(wmax1000_5a,5), 'k', linewidth=2)
# ax[2].plot(time/60, movmean(wmax1000_1b,5), 'b', linewidth=2)
# ax[2].plot(time/60, movmean(wmax1000_1c,5), 'r', linewidth=2)
# ax[2].set_title("Domain-maximum w(1 km)", fontsize=14)
# ax[2].set_xlabel('Time (min)', fontsize=12)
# ax[2].set_ylabel("w (m/s)", fontsize=12)
# ax[2].set_xlim([0,480])
# ax[2].set_ylim([0,25])

# for n in range(len(ax)):
#     ax[n].grid(visible=True, which='major', color='gray', linestyle='-')
#     ax[n].grid(visible=True, which='minor', color='lightgray', linestyle='-')
#     ax[n].xaxis.set_major_locator(MultipleLocator(30))
#     ax[n].xaxis.set_minor_locator(MultipleLocator(10))
#     ax[n].yaxis.set_major_locator(MultipleLocator(5))
#     ax[n].yaxis.set_minor_locator(MultipleLocator(1))
    
# if figsave:
#     plt.savefig(fp+'wmax_timeseries_ERA5_compare_3.png', dpi=300)
    
    


#%% Load translated swaths



fp = 'D:/brooks/era5-1_125m_final/'





ds = nc.Dataset(fp+'cm1out_000005.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
umove1 = ds.variables['umove'][:].data[0]
vmove1 = ds.variables['vmove'][:].data[0]
sws1 = ds.variables['sws2'][:].data[0,:,:]
shs1 = ds.variables['shs2'][:].data[0,:,:]
hail1 = ds.variables['hail2'][:].data[0,:,:]
dbz1 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000009.nc')
umove2 = ds.variables['umove'][:].data[0]
vmove2 = ds.variables['vmove'][:].data[0]
sws2 = ds.variables['sws2'][:].data[0,:,:]
shs2 = ds.variables['shs2'][:].data[0,:,:]
hail2 = ds.variables['hail2'][:].data[0,:,:]
dbz2 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000013.nc')
umove3 = ds.variables['umove'][:].data[0]
vmove3 = ds.variables['vmove'][:].data[0]
sws3 = ds.variables['sws2'][:].data[0,:,:]
shs3 = ds.variables['shs2'][:].data[0,:,:]
hail3 = ds.variables['hail2'][:].data[0,:,:]
dbz3 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000017.nc')
umove4 = ds.variables['umove'][:].data[0]
vmove4 = ds.variables['vmove'][:].data[0]
sws4 = ds.variables['sws2'][:].data[0,:,:]
shs4 = ds.variables['shs2'][:].data[0,:,:]
hail4 = ds.variables['hail2'][:].data[0,:,:]
dbz4 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000021.nc')
umove5 = ds.variables['umove'][:].data[0]
vmove5 = ds.variables['vmove'][:].data[0]
sws5 = ds.variables['sws2'][:].data[0,:,:]
shs5 = ds.variables['shs2'][:].data[0,:,:]
hail5 = ds.variables['hail2'][:].data[0,:,:]
dbz5 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000025.nc')
umove6 = ds.variables['umove'][:].data[0]
vmove6 = ds.variables['vmove'][:].data[0]
sws6 = ds.variables['sws2'][:].data[0,:,:]
shs6 = ds.variables['shs2'][:].data[0,:,:]
hail6 = ds.variables['hail2'][:].data[0,:,:]
dbz6 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000029.nc')
umove7 = ds.variables['umove'][:].data[0]
vmove7 = ds.variables['vmove'][:].data[0]
sws7 = ds.variables['sws2'][:].data[0,:,:]
shs7 = ds.variables['shs2'][:].data[0,:,:]
hail7 = ds.variables['hail2'][:].data[0,:,:]
dbz7 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000033.nc')
umove8 = ds.variables['umove'][:].data[0]
vmove8 = ds.variables['vmove'][:].data[0]
sws8 = ds.variables['sws2'][:].data[0,:,:]
shs8 = ds.variables['shs2'][:].data[0,:,:]
hail8 = ds.variables['hail2'][:].data[0,:,:]
dbz8 = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000037.nc')
umove9 = ds.variables['umove'][:].data[0]
vmove9 = ds.variables['vmove'][:].data[0]
sws9 = ds.variables['sws2'][:].data[0,:,:]
shs9 = ds.variables['shs2'][:].data[0,:,:]
hail9 = ds.variables['hail2'][:].data[0,:,:]
dbz9 = ds.variables['dbz'][:].data[0,0,:,:]
# wsp9 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove9, axis=0))**2 + 
#                (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove9, axis=0))**2)
ds.close()

ds = nc.Dataset(fp+'cm1out_000041.nc')
umove10 = ds.variables['umove'][:].data[0]
vmove10 = ds.variables['vmove'][:].data[0]
sws10 = ds.variables['sws2'][:].data[0,:,:]
shs10 = ds.variables['shs2'][:].data[0,:,:]
hail10 = ds.variables['hail2'][:].data[0,:,:]
dbz10 = ds.variables['dbz'][:].data[0,0,:,:]
# wsp10 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove10, axis=0))**2 + 
#                (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove10, axis=0))**2)
ds.close()

ds = nc.Dataset(fp+'cm1out_000045.nc')
umove11 = ds.variables['umove'][:].data[0]
vmove11 = ds.variables['vmove'][:].data[0]
sws11 = ds.variables['sws2'][:].data[0,:,:]
shs11 = ds.variables['shs2'][:].data[0,:,:]
hail11 = ds.variables['hail2'][:].data[0,:,:]
dbz11 = ds.variables['dbz'][:].data[0,0,:,:]
# wsp11 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove11, axis=0))**2 + 
#                (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove11, axis=0))**2)
ds.close()

ds = nc.Dataset(fp+'cm1out_000049.nc')
umove12 = ds.variables['umove'][:].data[0]
vmove12 = ds.variables['vmove'][:].data[0]
sws12 = ds.variables['sws2'][:].data[0,:,:]
shs12 = ds.variables['shs2'][:].data[0,:,:]
hail12 = ds.variables['hail2'][:].data[0,:,:]
dbz12 = ds.variables['dbz'][:].data[0,0,:,:]
# wsp12 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove12, axis=0))**2 + 
#                (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove12, axis=0))**2)
ds.close()


# x_added = umove1*3600/1000
# y_added = vmove1*3600/1000

# xh1 = xh + 100
# xh2 = xh1 + x_added
# xh3 = xh1 + 2*x_added
# xh4 = xh1 + 3*x_added
# xh5 = xh1 + 4*x_added
# xh6 = xh1 + 5*x_added
# xh7 = xh1 + 6*x_added
# xh8 = xh1 + 7*x_added
# xh9 = xh1 + 8*x_added
# xh10 = xh1 + 9*x_added
# # xh11 = xh1 + 10*x_added
# # xh12 = xh1 + 11*x_added

# yh1 = yh + 100
# yh2 = yh1 + y_added
# yh3 = yh1 + 2*y_added
# yh4 = yh1 + 3*y_added
# yh5 = yh1 + 4*y_added
# yh6 = yh1 + 5*y_added
# yh7 = yh1 + 6*y_added
# yh8 = yh1 + 7*y_added
# yh9 = yh1 + 8*y_added
# yh10 = yh1 + 9*y_added
# # yh11 = yh1 + 10*y_added
# # yh12 = yh1 + 11*y_added



xh1 = xh + 100
xh2 = xh1 + umove2*3600/1000
xh3 = xh2 + umove3*3600/1000
xh4 = xh3 + umove4*3600/1000
xh5 = xh4 + umove5*3600/1000
xh6 = xh5 + umove6*3600/1000
xh7 = xh6 + umove7*3600/1000
xh8 = xh7 + umove8*3600/1000
xh9 = xh8 + umove9*3600/1000
xh10 = xh9 + umove10*3600/1000
xh11 = xh10 + umove11*3600/1000
xh12 = xh11 + umove12*3600/1000

yh1 = yh + 100
yh2 = yh1 + vmove2*3600/1000
yh3 = yh2 + vmove3*3600/1000
yh4 = yh3 + vmove4*3600/1000
yh5 = yh4 + vmove5*3600/1000
yh6 = yh5 + vmove6*3600/1000
yh7 = yh6 + vmove7*3600/1000
yh8 = yh7 + vmove8*3600/1000
yh9 = yh8 + vmove9*3600/1000
yh10 = yh9 + vmove10*3600/1000
yh11 = yh10 + vmove11*3600/1000
yh12 = yh11 + vmove12*3600/1000

#%% Plot translated swaths


levs = [300,500]
cols = ['dimgray','k']
lws = [0.6,1]
dbz_levs = np.linspace(0,70,15)
# dbz_cols = ['limegreen','gold','darkorange','r']
dbz_lws = [1.25,1.25,1,1]


xt = [78, 120, 160, 202, 250, 295, 335, 380]
yt = [65, 73, 77, 78, 78, 78, 78, 78]


if 'era5-1_125m_test1' in fp:
    tstr = 'ERA5-1_test1'
    xt = [77, 120, 165, 215, 265, 320, 380, 435]
    yt = [73, 78, 80, 83, 89, 96, 102, 104]
elif 'era5-1_125m_test2' in fp:
    tstr = 'ERA5-1_test2'
    xt = [77, 120, 165, 215, 265, 325, 378, 425]
    yt = [70, 74, 77, 82, 87, 93, 94, 96]
elif 'era5-1_125m_test3' in fp:
    tstr = 'ERA5-1_test3'
    xt = [77, 116, 160, 210, 255, 306, 355, 405]
    yt = [71, 75, 79, 79, 80, 81, 82, 86]
elif 'era5-1_125m_test4' in fp:
    tstr = 'ERA5-1_test4'
    xt = [75, 115, 157, 201, 247, 298, 348, 397]
    yt = [74, 83, 85, 85, 85, 85, 85, 85]
elif 'era5-1_125m_test5' in fp:
    tstr = 'ERA5-1_test5'
    xt = [78, 120, 160, 202, 250, 295, 335, 380]
    yt = [65, 73, 73, 75, 75, 78, 79, 79]
elif 'era5-1_125m_test6' in fp:
    tstr = 'ERA5-1_test6'
    xt = [80, 122, 163, 205, 250, 292, 338, 380]
    yt = [60, 67, 68, 70, 71, 73, 73, 73]
elif 'era5-1_125m_test7' in fp:
    tstr = 'ERA5-1_test7'
    xt = [78, 120, 165, 208, 255, 300, 342, 383]
    yt = [65, 68, 72, 73, 75, 75, 73, 72]
elif 'era5-1_125m_test8' in fp:
    tstr = 'ERA5-1_test8'
    xt = [78, 120, 160, 202, 
          245, 290, 335, 375, 
          420, 460, 510, 565]
    yt = [67, 72, 75, 78, 
          78, 79, 79, 80, 
          80, 81, 82, 83]
elif 'era5-1_125m_final' in fp:
    tstr = 'ERA5-1_final'
    xt = [78, 120, 160, 202,
          250, 295, 335, 385,
          435, 480, 530, 580]
    yt = [65, 65, 65, 65,
          65, 65, 65, 65,
          65, 65, 65, 65]
    
    
# xl = [50,500]
yl = [50,200]
xl = [50,675]

if 'era5-1' in fp:
    sounding_str = 'ERA5 2025-08-20-T21z (50.75,-114.0)'
elif 'era5-2' in fp:
    sounding_str = 'ERA5 2025-08-20-T21z (51.0,-114.25)'
elif 'era5-3' in fp:
    sounding_str = 'ERA5 2025-08-20-T21z (51.25,-113.25)'
elif 'hrdps' in fp:
    sounding_str = 'HRDPS 2025-08-20-T21z (51.166,-113.135)'



figsave = False



# fig,ax = plt.subplots(1, 1, figsize=(8.5,2.75), subplot_kw=dict(aspect=1), layout='constrained')
fig,ax = plt.subplots(1, 1, figsize=(10,2.5), subplot_kw=dict(aspect=1), layout='constrained')

c = ax.contourf(xh1, yh1, np.ma.masked_array(dbz1, dbz1<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh2, yh2, np.ma.masked_array(dbz2, dbz2<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh3, yh3, np.ma.masked_array(dbz3, dbz3<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh4, yh4, np.ma.masked_array(dbz4, dbz4<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh5, yh5, np.ma.masked_array(dbz5, dbz5<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh6, yh6, np.ma.masked_array(dbz6, dbz6<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh7, yh7, np.ma.masked_array(dbz7, dbz7<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh8, yh8, np.ma.masked_array(dbz8, dbz8<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh9, yh9, np.ma.masked_array(dbz9, dbz9<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh10, yh10, np.ma.masked_array(dbz10, dbz10<20), levels=dbz_levs, vmin=0, vmax=70)
# ax.contourf(xh11, yh11, np.ma.masked_array(dbz11, dbz11<20), levels=dbz_levs, vmin=0, vmax=70)
# ax.contourf(xh12, yh12, np.ma.masked_array(dbz12, dbz12<20), levels=dbz_levs, vmin=0, vmax=70)

cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_ticks(np.linspace(0,70,8))
cb.set_label('Radar reflectivity (dBZ)', fontsize=10)

ax.contour(xh1, yh1, shs1, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh2, yh2, shs2, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh3, yh3, shs3, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh4, yh4, shs4, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh5, yh5, shs5, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh6, yh6, shs6, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh7, yh7, shs7, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh8, yh8, shs8, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh9, yh9, shs9, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh10, yh10, shs10, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh11, yh11, shs11, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh12, yh12, shs12, levels=levs, colors=cols, linewidths=lws)
ax.set_xlim(xl)
ax.set_ylim(yl)
ax.set_xlabel('Translated x (km)', fontsize=10)
ax.set_ylabel('y (km)', fontsize=10)
ax.set_title(f"{sounding_str} - Updraft helicity swaths (0-12 h)", fontsize=10)
l1, = ax.plot([-2,-1], [-2,-1], 'gray', linewidth=0.75)
l2, = ax.plot([-2,-1], [-2,-1], 'k', linewidth=1)
ax.legend(handles=[l1,l2], labels=['300 m2/s2','500 m2/s2'], loc='lower right', fontsize=8)

ax.text(xt[0], yt[0], '1 h', fontsize=9, fontweight='bold')
ax.text(xt[1], yt[1], '2 h', fontsize=9, fontweight='bold')
ax.text(xt[2], yt[2], '3 h', fontsize=9, fontweight='bold')
ax.text(xt[3], yt[3], '4 h', fontsize=9, fontweight='bold')
ax.text(xt[4], yt[4], '5 h', fontsize=9, fontweight='bold')
ax.text(xt[5], yt[5], '6 h', fontsize=9, fontweight='bold')
ax.text(xt[6], yt[6], '7 h', fontsize=9, fontweight='bold')
ax.text(xt[7], yt[7], '8 h', fontsize=9, fontweight='bold')
ax.text(xt[8], yt[8], '9 h', fontsize=9, fontweight='bold')
ax.text(xt[9], yt[9], '10 h', fontsize=9, fontweight='bold')
# ax.text(xt[10], yt[10], '11 h', fontsize=9, fontweight='bold')
# ax.text(xt[11], yt[11], '12 h', fontsize=9, fontweight='bold')

if figsave:
    plt.savefig(fp+f"figs/dbz_uh_swath.png", dpi=300)






levs = [0.01,0.1]
cols = ['gray','k']
lws = [0.75,1]

# fig,ax = plt.subplots(1, 1, figsize=(8.5,2.75), subplot_kw=dict(aspect=1), layout='constrained')
fig,ax = plt.subplots(1, 1, figsize=(10,2.5), subplot_kw=dict(aspect=1), layout='constrained')

c = ax.contourf(xh1, yh1, np.ma.masked_array(dbz1, dbz1<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh2, yh2, np.ma.masked_array(dbz2, dbz2<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh3, yh3, np.ma.masked_array(dbz3, dbz3<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh4, yh4, np.ma.masked_array(dbz4, dbz4<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh5, yh5, np.ma.masked_array(dbz5, dbz5<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh6, yh6, np.ma.masked_array(dbz6, dbz6<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh7, yh7, np.ma.masked_array(dbz7, dbz7<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh8, yh8, np.ma.masked_array(dbz8, dbz8<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh9, yh9, np.ma.masked_array(dbz9, dbz9<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh10, yh10, np.ma.masked_array(dbz10, dbz10<20), levels=dbz_levs, vmin=0, vmax=70)
# ax.contourf(xh11, yh11, np.ma.masked_array(dbz11, dbz11<20), levels=dbz_levs, vmin=0, vmax=70)
# ax.contourf(xh12, yh12, np.ma.masked_array(dbz12, dbz12<20), levels=dbz_levs, vmin=0, vmax=70)

cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_ticks(np.linspace(0,70,8))
cb.set_label('Radar reflectivity (dBZ)', fontsize=10)

ax.contour(xh1, yh1, hail1, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh2, yh2, hail2, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh3, yh3, hail3, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh4, yh4, hail4, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh5, yh5, hail5, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh6, yh6, hail6, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh7, yh7, hail7, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh8, yh8, hail8, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh9, yh9, hail9, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh10, yh10, hail10, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh11, yh11, hail11, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh12, yh12, hail12, levels=levs, colors=cols, linewidths=lws)
ax.set_xlim(xl)
ax.set_ylim(yl)
ax.set_xlabel('Translated x (km)', fontsize=10)
ax.set_ylabel('y (km)', fontsize=10)
ax.set_title(f"{sounding_str} - Accumulated hail swaths (0-12 h)", fontsize=10)
l1, = ax.plot([-2,-1], [-2,-1], 'gray', linewidth=0.75)
l2, = ax.plot([-2,-1], [-2,-1], 'k', linewidth=1)
ax.legend(handles=[l1,l2], labels=['0.1 mm','1 mm'], loc='lower right', fontsize=8)

ax.text(xt[0], yt[0], '1 h', fontsize=9, fontweight='bold')
ax.text(xt[1], yt[1], '2 h', fontsize=9, fontweight='bold')
ax.text(xt[2], yt[2], '3 h', fontsize=9, fontweight='bold')
ax.text(xt[3], yt[3], '4 h', fontsize=9, fontweight='bold')
ax.text(xt[4], yt[4], '5 h', fontsize=9, fontweight='bold')
ax.text(xt[5], yt[5], '6 h', fontsize=9, fontweight='bold')
ax.text(xt[6], yt[6], '7 h', fontsize=9, fontweight='bold')
ax.text(xt[7], yt[7], '8 h', fontsize=9, fontweight='bold')
ax.text(xt[8], yt[8], '9 h', fontsize=9, fontweight='bold')
ax.text(xt[9], yt[9], '10 h', fontsize=9, fontweight='bold')
# ax.text(xt[10], yt[10], '11 h', fontsize=9, fontweight='bold')
# ax.text(xt[11], yt[11], '12 h', fontsize=9, fontweight='bold')

if figsave:
    plt.savefig(fp+f"figs/dbz_hail_swath.png", dpi=300)






levs = [20,25.7]
cols = ['gray','k']
lws = [0.75,1]

# fig,ax = plt.subplots(1, 1, figsize=(8.5,2.75), subplot_kw=dict(aspect=1), layout='constrained')
fig,ax = plt.subplots(1, 1, figsize=(10,2.5), subplot_kw=dict(aspect=1), layout='constrained')

c = ax.contourf(xh1, yh1, np.ma.masked_array(dbz1, dbz1<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh2, yh2, np.ma.masked_array(dbz2, dbz2<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh3, yh3, np.ma.masked_array(dbz3, dbz3<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh4, yh4, np.ma.masked_array(dbz4, dbz4<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh5, yh5, np.ma.masked_array(dbz5, dbz5<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh6, yh6, np.ma.masked_array(dbz6, dbz6<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh7, yh7, np.ma.masked_array(dbz7, dbz7<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh8, yh8, np.ma.masked_array(dbz8, dbz8<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh9, yh9, np.ma.masked_array(dbz9, dbz9<20), levels=dbz_levs, vmin=0, vmax=70)
ax.contourf(xh10, yh10, np.ma.masked_array(dbz10, dbz10<20), levels=dbz_levs, vmin=0, vmax=70)
# ax.contourf(xh11, yh11, np.ma.masked_array(dbz11, dbz11<20), levels=dbz_levs, vmin=0, vmax=70)
# ax.contourf(xh12, yh12, np.ma.masked_array(dbz12, dbz12<20), levels=dbz_levs, vmin=0, vmax=70)

cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_ticks(np.linspace(0,70,8))
cb.set_label('Radar reflectivity (dBZ)', fontsize=10)

ax.contour(xh1, yh1, sws1, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh2, yh2, sws2, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh3, yh3, sws3, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh4, yh4, sws4, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh5, yh5, sws5, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh6, yh6, sws6, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh7, yh7, sws7, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh8, yh8, sws8, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh9, yh9, sws9, levels=levs, colors=cols, linewidths=lws)
ax.contour(xh10, yh10, sws10, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh11, yh11, sws11, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh12, yh12, sws12, levels=levs, colors=cols, linewidths=lws)
ax.set_xlim(xl)
ax.set_ylim(yl)
ax.set_xlabel('Translated x (km)', fontsize=10)
ax.set_ylabel('y (km)', fontsize=10)
ax.set_title(f"{sounding_str} - Surface wind swaths (0-12 h)", fontsize=10)
l1, = ax.plot([-2,-1], [-2,-1], 'gray', linewidth=0.75)
l2, = ax.plot([-2,-1], [-2,-1], 'k', linewidth=1)
ax.legend(handles=[l1,l2], labels=['20 m/s','26 m/s'], loc='lower right', fontsize=8)

ax.text(xt[0], yt[0], '1 h', fontsize=9, fontweight='bold')
ax.text(xt[1], yt[1], '2 h', fontsize=9, fontweight='bold')
ax.text(xt[2], yt[2], '3 h', fontsize=9, fontweight='bold')
ax.text(xt[3], yt[3], '4 h', fontsize=9, fontweight='bold')
ax.text(xt[4], yt[4], '5 h', fontsize=9, fontweight='bold')
ax.text(xt[5], yt[5], '6 h', fontsize=9, fontweight='bold')
ax.text(xt[6], yt[6], '7 h', fontsize=9, fontweight='bold')
ax.text(xt[7], yt[7], '8 h', fontsize=9, fontweight='bold')
ax.text(xt[8], yt[8], '9 h', fontsize=9, fontweight='bold')
ax.text(xt[9], yt[9], '10 h', fontsize=9, fontweight='bold')
# ax.text(xt[10], yt[10], '11 h', fontsize=9, fontweight='bold')
# ax.text(xt[11], yt[11], '12 h', fontsize=9, fontweight='bold')

if figsave:
    plt.savefig(fp+f"figs/dbz_wind_swath.png", dpi=300)



#%% Swath plots but just 6-12 h -- Load data

fp = 'D:/brooks/era5-1_125m_final/'



ds = nc.Dataset(fp+'cm1out_000025.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['zh'][:].data
iz80 = np.argmin(abs(zh-0.07))

umove1 = ds.variables['umove'][:].data[0]
vmove1 = ds.variables['vmove'][:].data[0]
sws1 = ds.variables['sws2'][:].data[0,:,:]
shs1 = ds.variables['shs2'][:].data[0,:,:]
hail1 = ds.variables['hail2'][:].data[0,:,:]
dbz1 = ds.variables['dbz'][:].data[0,0,:,:]
# wsp1 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove1, axis=0))**2 + 
#                (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove1, axis=0))**2)
# dmi1 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,:,:,:], ds.variables['p3_dmi2'][:].data[0,:,:,:],
#                           ds.variables['p3_dmi3'][:].data[0,:,:,:], ds.variables['p3_dmi4'][:].data[0,:,:,:]]), axis=0)
# sfcdmi1 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,0,:,:], ds.variables['p3_dmi2'][:].data[0,0,:,:],
#                              ds.variables['p3_dmi3'][:].data[0,0,:,:], ds.variables['p3_dmi4'][:].data[0,0,:,:]]), axis=0)
ds.close()

ds = nc.Dataset(fp+'cm1out_000029.nc')
umove2 = ds.variables['umove'][:].data[0]
vmove2 = ds.variables['vmove'][:].data[0]
sws2 = ds.variables['sws2'][:].data[0,:,:]
shs2 = ds.variables['shs2'][:].data[0,:,:]
hail2 = ds.variables['hail2'][:].data[0,:,:]
dbz2 = ds.variables['dbz'][:].data[0,0,:,:]
# wsp2 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove2, axis=0))**2 + 
#                (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove2, axis=0))**2)
# dmi2 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,:,:,:], ds.variables['p3_dmi2'][:].data[0,:,:,:],
#                           ds.variables['p3_dmi3'][:].data[0,:,:,:], ds.variables['p3_dmi4'][:].data[0,:,:,:]]), axis=0)
# sfcdmi2 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,0,:,:], ds.variables['p3_dmi2'][:].data[0,0,:,:],
#                              ds.variables['p3_dmi3'][:].data[0,0,:,:], ds.variables['p3_dmi4'][:].data[0,0,:,:]]), axis=0)
ds.close()

ds = nc.Dataset(fp+'cm1out_000033.nc')
umove3 = ds.variables['umove'][:].data[0]
vmove3 = ds.variables['vmove'][:].data[0]
sws3 = ds.variables['sws2'][:].data[0,:,:]
shs3 = ds.variables['shs2'][:].data[0,:,:]
hail3 = ds.variables['hail2'][:].data[0,:,:]
dbz3 = ds.variables['dbz'][:].data[0,0,:,:]
# wsp3 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove3, axis=0))**2 + 
#                (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove3, axis=0))**2)
# dmi3 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,:,:,:], ds.variables['p3_dmi2'][:].data[0,:,:,:],
#                           ds.variables['p3_dmi3'][:].data[0,:,:,:], ds.variables['p3_dmi4'][:].data[0,:,:,:]]), axis=0)
# sfcdmi3 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,0,:,:], ds.variables['p3_dmi2'][:].data[0,0,:,:],
#                              ds.variables['p3_dmi3'][:].data[0,0,:,:], ds.variables['p3_dmi4'][:].data[0,0,:,:]]), axis=0)
ds.close()

ds = nc.Dataset(fp+'cm1out_000037.nc')
umove4 = ds.variables['umove'][:].data[0]
vmove4 = ds.variables['vmove'][:].data[0]
sws4 = ds.variables['sws2'][:].data[0,:,:]
shs4 = ds.variables['shs2'][:].data[0,:,:]
hail4 = ds.variables['hail2'][:].data[0,:,:]
dbz4 = ds.variables['dbz'][:].data[0,0,:,:]
# wsp4 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove4, axis=0))**2 + 
#                (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove4, axis=0))**2)
# dmi4 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,:,:,:], ds.variables['p3_dmi2'][:].data[0,:,:,:],
#                           ds.variables['p3_dmi3'][:].data[0,:,:,:], ds.variables['p3_dmi4'][:].data[0,:,:,:]]), axis=0)
# sfcdmi4 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,0,:,:], ds.variables['p3_dmi2'][:].data[0,0,:,:],
#                              ds.variables['p3_dmi3'][:].data[0,0,:,:], ds.variables['p3_dmi4'][:].data[0,0,:,:]]), axis=0)
ds.close()

ds = nc.Dataset(fp+'cm1out_000041.nc')
umove5 = ds.variables['umove'][:].data[0]
vmove5 = ds.variables['vmove'][:].data[0]
sws5 = ds.variables['sws2'][:].data[0,:,:]
shs5 = ds.variables['shs2'][:].data[0,:,:]
hail5 = ds.variables['hail2'][:].data[0,:,:]
dbz5 = ds.variables['dbz'][:].data[0,0,:,:]
# wsp5 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove5, axis=0))**2 + 
#                (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove5, axis=0))**2)
# dmi5 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,:,:,:], ds.variables['p3_dmi2'][:].data[0,:,:,:],
#                           ds.variables['p3_dmi3'][:].data[0,:,:,:], ds.variables['p3_dmi4'][:].data[0,:,:,:]]), axis=0)
# sfcdmi5 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,0,:,:], ds.variables['p3_dmi2'][:].data[0,0,:,:],
#                              ds.variables['p3_dmi3'][:].data[0,0,:,:], ds.variables['p3_dmi4'][:].data[0,0,:,:]]), axis=0)
ds.close()

ds = nc.Dataset(fp+'cm1out_000045.nc')
umove6 = ds.variables['umove'][:].data[0]
vmove6 = ds.variables['vmove'][:].data[0]
sws6 = ds.variables['sws2'][:].data[0,:,:]
shs6 = ds.variables['shs2'][:].data[0,:,:]
hail6 = ds.variables['hail2'][:].data[0,:,:]
dbz6 = ds.variables['dbz'][:].data[0,0,:,:]
# wsp6 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove6, axis=0))**2 + 
#                (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove6, axis=0))**2)
# dmi6 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,:,:,:], ds.variables['p3_dmi2'][:].data[0,:,:,:],
#                           ds.variables['p3_dmi3'][:].data[0,:,:,:], ds.variables['p3_dmi4'][:].data[0,:,:,:]]), axis=0)
# sfcdmi6 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,0,:,:], ds.variables['p3_dmi2'][:].data[0,0,:,:],
#                              ds.variables['p3_dmi3'][:].data[0,0,:,:], ds.variables['p3_dmi4'][:].data[0,0,:,:]]), axis=0)
ds.close()

ds = nc.Dataset(fp+'cm1out_000049.nc')
umove7 = ds.variables['umove'][:].data[0]
vmove7 = ds.variables['vmove'][:].data[0]
sws7 = ds.variables['sws2'][:].data[0,:,:]
shs7 = ds.variables['shs2'][:].data[0,:,:]
hail7 = ds.variables['hail2'][:].data[0,:,:]
dbz7 = ds.variables['dbz'][:].data[0,0,:,:]
# wsp7 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove7, axis=0))**2 + 
#                (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove7, axis=0))**2)
# dmi7 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,:,:,:], ds.variables['p3_dmi2'][:].data[0,:,:,:],
#                           ds.variables['p3_dmi3'][:].data[0,:,:,:], ds.variables['p3_dmi4'][:].data[0,:,:,:]]), axis=0)
# sfcdmi7 = np.max(np.asarray([ds.variables['p3_dmi1'][:].data[0,0,:,:], ds.variables['p3_dmi2'][:].data[0,0,:,:],
#                              ds.variables['p3_dmi3'][:].data[0,0,:,:], ds.variables['p3_dmi4'][:].data[0,0,:,:]]), axis=0)
ds.close()





dbfile = open(fp+'pkls/wind_mechanisms_360min.pkl', 'rb'); crit1 = pickle.load(dbfile); dbfile.close()
dbfile = open(fp+'pkls/wind_mechanisms_420min.pkl', 'rb'); crit2 = pickle.load(dbfile); dbfile.close()
dbfile = open(fp+'pkls/wind_mechanisms_480min.pkl', 'rb'); crit3 = pickle.load(dbfile); dbfile.close()
dbfile = open(fp+'pkls/wind_mechanisms_540min.pkl', 'rb'); crit4 = pickle.load(dbfile); dbfile.close()
dbfile = open(fp+'pkls/wind_mechanisms_600min.pkl', 'rb'); crit5 = pickle.load(dbfile); dbfile.close()
dbfile = open(fp+'pkls/wind_mechanisms_660min.pkl', 'rb'); crit6 = pickle.load(dbfile); dbfile.close()
dbfile = open(fp+'pkls/wind_mechanisms_720min.pkl', 'rb'); crit7 = pickle.load(dbfile); dbfile.close()


crit = {'1':crit1, '2':crit2, '3':crit3, '4':crit4, '5':crit5, '6':crit6, '7':crit7}




xh1 = xh + 50
xh2 = xh1 + umove2*3600/1000
xh3 = xh2 + umove3*3600/1000
xh4 = xh3 + umove4*3600/1000
xh5 = xh4 + umove5*3600/1000
xh6 = xh5 + umove6*3600/1000
xh7 = xh6 + umove7*3600/1000

yh1 = yh + 50
yh2 = yh1 + vmove2*3600/1000
yh3 = yh2 + vmove3*3600/1000
yh4 = yh3 + vmove4*3600/1000
yh5 = yh4 + vmove5*3600/1000
yh6 = yh5 + vmove6*3600/1000
yh7 = yh6 + vmove7*3600/1000

xx = {'1':xh1, '2':xh2, '3':xh3, '4':xh4, '5':xh5, '6':xh6, '7':xh7}
yy = {'1':yh1, '2':yh2, '3':yh3, '4':yh4, '5':yh5, '6':yh6, '7':yh7}


#%% Swath plots but just 6-12 h -- Make plots

dbz_levs = np.linspace(0,70,15); dbz_cm = "HomeyerRainbow"
wsp_levs = np.linspace(0,35,36); wsp_cm = "Blues"
shs_levs = [500]; shs_cols = ['k']; shs_lws = [0.6]
hail_levs = [0.1]; hail_cols = ['k']; hail_lws = [0.6]
sws_levs = [25.7]; sws_cols = ['k']; sws_lws = [0.6]


V_thres = 25.7

tstr = 'ERA5-1_final'
xt = [46, 86, 136, 186, 231, 281, 331]
yt = [5, 5, 5, 5, 5, 5, 5]


xl = [0,400]
yl = [0,150]



figsave = False



### dbz + UH
fig,ax = plt.subplots(1, 1, figsize=(7.5,2.75), subplot_kw=dict(aspect=1), layout='constrained')

c = ax.contourf(xh1, yh1, np.ma.masked_array(dbz1, dbz1<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh2, yh2, np.ma.masked_array(dbz2, dbz2<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh3, yh3, np.ma.masked_array(dbz3, dbz3<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh4, yh4, np.ma.masked_array(dbz4, dbz4<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh5, yh5, np.ma.masked_array(dbz5, dbz5<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh6, yh6, np.ma.masked_array(dbz6, dbz6<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh7, yh7, np.ma.masked_array(dbz7, dbz7<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)

cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_ticks(np.linspace(0,70,8))
cb.set_label('Radar reflectivity (dBZ)', fontsize=10)

ax.contour(xh1[:-50], yh1, shs1[:,:-50], levels=shs_levs, colors=shs_cols, linewidths=shs_lws)
ax.contour(xh2[:-50], yh2, shs2[:,:-50], levels=shs_levs, colors=shs_cols, linewidths=shs_lws)
ax.contour(xh3[:-50], yh3, shs3[:,:-50], levels=shs_levs, colors=shs_cols, linewidths=shs_lws)
ax.contour(xh4[:-50], yh4, shs4[:,:-50], levels=shs_levs, colors=shs_cols, linewidths=shs_lws)
ax.contour(xh5[:-50], yh5, shs5[:,:-50], levels=shs_levs, colors=shs_cols, linewidths=shs_lws)
ax.contour(xh6[:-50], yh6, shs6[:,:-50], levels=shs_levs, colors=shs_cols, linewidths=shs_lws)
ax.contour(xh7[:-50], yh7, shs7[:,:-50], levels=shs_levs, colors=shs_cols, linewidths=shs_lws)

ax.set_xlim(xl)
ax.set_ylim(yl)
ax.set_xlabel('Translated x (km)', fontsize=10)
ax.set_ylabel('Translated y (km)', fontsize=10)
ax.set_title(f"Surface reflectivity, updraft helicity swaths (6-12 h)", fontsize=10)
l1, = ax.plot([-2,-1], [-2,-1], 'dimgray', linewidth=1)
l2, = ax.plot([-2,-1], [-2,-1], 'k', linewidth=1)
ax.legend(handles=[l2], labels=['UH=500 m2/s2'], loc='upper right', fontsize=8)

ax.text(xt[0], yt[0], '6 h', fontsize=9, fontweight='bold')
ax.text(xt[1], yt[1], '7 h', fontsize=9, fontweight='bold')
ax.text(xt[2], yt[2], '8 h', fontsize=9, fontweight='bold')
ax.text(xt[3], yt[3], '9 h', fontsize=9, fontweight='bold')
ax.text(xt[4], yt[4], '10 h', fontsize=9, fontweight='bold')
ax.text(xt[5], yt[5], '11 h', fontsize=9, fontweight='bold')
ax.text(xt[6], yt[6], '12 h', fontsize=9, fontweight='bold')

if figsave:
    plt.savefig(fp+f"figs/dbz_uh_swath_6-12H.png", dpi=300)




### dbz + hail
fig,ax = plt.subplots(1, 1, figsize=(7.5,2.75), subplot_kw=dict(aspect=1), layout='constrained')

c = ax.contourf(xh1, yh1, np.ma.masked_array(dbz1, dbz1<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh2, yh2, np.ma.masked_array(dbz2, dbz2<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh3, yh3, np.ma.masked_array(dbz3, dbz3<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh4, yh4, np.ma.masked_array(dbz4, dbz4<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh5, yh5, np.ma.masked_array(dbz5, dbz5<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh6, yh6, np.ma.masked_array(dbz6, dbz6<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh7, yh7, np.ma.masked_array(dbz7, dbz7<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)

cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_ticks(np.linspace(0,70,8))
cb.set_label('Radar reflectivity (dBZ)', fontsize=10)

ax.contour(xh1, yh1, hail1, levels=hail_levs, colors=hail_cols, linewidths=hail_lws)
ax.contour(xh2, yh2, hail2, levels=hail_levs, colors=hail_cols, linewidths=hail_lws)
ax.contour(xh3, yh3, hail3, levels=hail_levs, colors=hail_cols, linewidths=hail_lws)
ax.contour(xh4, yh4, hail4, levels=hail_levs, colors=hail_cols, linewidths=hail_lws)
ax.contour(xh5, yh5, hail5, levels=hail_levs, colors=hail_cols, linewidths=hail_lws)
ax.contour(xh6, yh6, hail6, levels=hail_levs, colors=hail_cols, linewidths=hail_lws)
ax.contour(xh7, yh7, hail7, levels=hail_levs, colors=hail_cols, linewidths=hail_lws)

ax.set_xlim(xl)
ax.set_ylim(yl)
ax.set_xlabel('Translated x (km)', fontsize=10)
ax.set_ylabel('Translated y (km)', fontsize=10)
ax.set_title(f"Surface reflectivity, accumulated surface hail swaths (6-12 h)", fontsize=10)
l1, = ax.plot([-2,-1], [-2,-1], 'dimgray', linewidth=1)
l2, = ax.plot([-2,-1], [-2,-1], 'k', linewidth=1)
ax.legend(handles=[l2], labels=['1 mm'], loc='upper right', fontsize=8)

ax.text(xt[0], yt[0], '6 h', fontsize=9, fontweight='bold')
ax.text(xt[1], yt[1], '7 h', fontsize=9, fontweight='bold')
ax.text(xt[2], yt[2], '8 h', fontsize=9, fontweight='bold')
ax.text(xt[3], yt[3], '9 h', fontsize=9, fontweight='bold')
ax.text(xt[4], yt[4], '10 h', fontsize=9, fontweight='bold')
ax.text(xt[5], yt[5], '11 h', fontsize=9, fontweight='bold')
ax.text(xt[6], yt[6], '12 h', fontsize=9, fontweight='bold')

if figsave:
    plt.savefig(fp+f"figs/dbz_hail_swath_6-12H.png", dpi=300)




### dbz + wpsd
fig,ax = plt.subplots(1, 1, figsize=(7.5,2.75), subplot_kw=dict(aspect=1), layout='constrained')

c = ax.contourf(xh1, yh1, np.ma.masked_array(dbz1, dbz1<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh2, yh2, np.ma.masked_array(dbz2, dbz2<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh3, yh3, np.ma.masked_array(dbz3, dbz3<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh4, yh4, np.ma.masked_array(dbz4, dbz4<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh5, yh5, np.ma.masked_array(dbz5, dbz5<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh6, yh6, np.ma.masked_array(dbz6, dbz6<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)
ax.contourf(xh7, yh7, np.ma.masked_array(dbz7, dbz7<20), levels=dbz_levs, vmin=0, vmax=70, cmap=dbz_cm)

cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_ticks(np.linspace(0,70,8))
cb.set_label('Radar reflectivity (dBZ)', fontsize=10)

ax.contour(xh1, yh1, sws1, levels=sws_levs, colors=sws_cols, linewidths=sws_lws)
ax.contour(xh2, yh2, sws2, levels=sws_levs, colors=sws_cols, linewidths=sws_lws)
ax.contour(xh3, yh3, sws3, levels=sws_levs, colors=sws_cols, linewidths=sws_lws)
ax.contour(xh4, yh4, sws4, levels=sws_levs, colors=sws_cols, linewidths=sws_lws)
ax.contour(xh5, yh5, sws5, levels=sws_levs, colors=sws_cols, linewidths=sws_lws)
ax.contour(xh6, yh6, sws6, levels=sws_levs, colors=sws_cols, linewidths=sws_lws)
ax.contour(xh7, yh7, sws7, levels=sws_levs, colors=sws_cols, linewidths=sws_lws)

ax.set_xlim(xl)
ax.set_ylim(yl)
ax.set_xlabel('Translated x (km)', fontsize=10)
ax.set_ylabel('Translated y (km)', fontsize=10)
ax.set_title(f"Surface reflectivity, surface wind swaths (6-12 h)", fontsize=10)
l1, = ax.plot([-2,-1], [-2,-1], 'dimgray', linewidth=1)
l2, = ax.plot([-2,-1], [-2,-1], 'k', linewidth=1)
ax.legend(handles=[l2], labels=['V=25.7 m/s'], loc='upper right', fontsize=8)

ax.text(xt[0], yt[0], '6 h', fontsize=9, fontweight='bold')
ax.text(xt[1], yt[1], '7 h', fontsize=9, fontweight='bold')
ax.text(xt[2], yt[2], '8 h', fontsize=9, fontweight='bold')
ax.text(xt[3], yt[3], '9 h', fontsize=9, fontweight='bold')
ax.text(xt[4], yt[4], '10 h', fontsize=9, fontweight='bold')
ax.text(xt[5], yt[5], '11 h', fontsize=9, fontweight='bold')
ax.text(xt[6], yt[6], '12 h', fontsize=9, fontweight='bold')

if figsave:
    plt.savefig(fp+f"figs/dbz_wind_swath_6-12H.png", dpi=300)




#%% Plot base state sounding

import xarray as xr
import metpy.calc as mc
from metpy.units import units
from metpy.plots import SkewT,Hodograph
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable




figsave = False



ds = xr.open_dataset("D:/brooks/rerun1/cm1out_000001.nc")
time = ds['time'][:].values[0]
xh = ds['xh'][:].values
yh = ds['yh'][:].values
zh = ds['zh'][:].values
umove = ds['umove'].values[0]
vmove = ds['vmove'].values[0]

u0 = ds.sel(time=time, yh=yh[0], xh=xh[0])['uinterp'].values + umove
v0 = ds.sel(time=time, yh=yh[0], xh=xh[0])['vinterp'].values + vmove
th0 = ds.sel(time=time, yh=yh[0], xh=xh[0])['th'].values
qv0 = ds.sel(time=time, yh=yh[0], xh=xh[0])['qv'].values
prs0 = ds.sel(time=time, yh=yh[0], xh=xh[0])['prs'].values/100
ds.close()


T0 = th0 * (prs0/1000)**0.286
e0 = (qv0*prs0) / (0.622+qv0)
Td0 = 243.5 / ((17.67 / (np.log(e0/6.112))) - 1) + 273.15
T0_parcel = mc.parcel_profile(prs0*units.hPa, T0[0]*units.K, Td0[0]*units.K)
lcl_pressure,lcl_temperature = mc.lcl(prs0[0]*units.hPa, T0[0]*units.K, Td0[0]*units.K)
ip = np.where(prs0<=lcl_pressure.magnitude)[0][0]
dp = prs0[ip-1] - prs0[ip]
w1 = (lcl_pressure.magnitude - prs0[ip])/dp
w2 = (prs0[ip-1] - lcl_pressure.magnitude)/dp
lcl_height = w1*zh[ip-1] + w2*zh[ip]

[rm,lm,mm] = mc.bunkers_storm_motion(prs0*units.hPa, u0*units('m/s'), v0*units('m/s'), zh*units.km)
mw = mc.wind_speed(mm[0], mm[1]).magnitude
s = mc.bulk_shear(prs0*units.hPa, u0*units('m/s'), v0*units('m/s'), zh*1000*units.m, depth=6000*units.m)
s06 = np.sqrt(s[0].magnitude**2 + s[1].magnitude**2)
s = mc.bulk_shear(prs0*units.hPa, u0*units('m/s'), v0*units('m/s'), zh*1000*units.m, depth=3000*units.m)
s03 = np.sqrt(s[0].magnitude**2 + s[1].magnitude**2)
s = mc.bulk_shear(prs0*units.hPa, u0*units('m/s'), v0*units('m/s'), zh*1000*units.m, depth=1000*units.m)
s01 = np.sqrt(s[0].magnitude**2 + s[1].magnitude**2)
s = mc.storm_relative_helicity(zh*1000*units.m, u0*units('m/s'), v0*units('m/s'), depth=1000*units.m, storm_u=rm[0], storm_v=rm[1])
srh01 = s[2].magnitude
s = mc.storm_relative_helicity(zh*1000*units.m, u0*units('m/s'), v0*units('m/s'), depth=3000*units.m, storm_u=rm[0], storm_v=rm[1])
srh03 = s[2].magnitude

mlcape,mlcin = mc.mixed_layer_cape_cin(prs0*units.hPa, T0*units.K, Td0*units.K)
mucape,mucin = mc.most_unstable_cape_cin(prs0*units.hPa, T0*units.K, Td0*units.K)
sbcape,sbcin = mc.surface_based_cape_cin(prs0*units.hPa, T0*units.K, Td0*units.K)


u_barbs = []
v_barbs = []
p_barbs = []
for p in np.arange(1000,50,-50):
    ip = np.where(np.abs(prs0-p) == np.min(np.abs(prs0-p)))[0][0]
    u_barbs.append(u0[ip])
    v_barbs.append(v0[ip])
    p_barbs.append(prs0[ip])
u_barbs = (np.array(u_barbs)*units('m/s')).to('knots')
v_barbs = (np.array(v_barbs)*units('m/s')).to('knots')
p_barbs = np.array(p_barbs)




fig = plt.figure(figsize=(10,9))
skew = SkewT(fig=fig)
skew.plot(prs0, (T0-273.15), '-r', linewidth=2)
skew.plot(prs0, (Td0-273.15), '-g', linewidth=2)
skew.plot(prs0, np.array(T0_parcel.magnitude[:])-273.15, '-k', linewidth=2)

skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black') # Plot LCL as black dot
skew.plot_dry_adiabats(linewidths=1)
skew.plot_moist_adiabats(linewidths=1)
skew.plot_mixing_lines()
b = skew.plot_barbs(p_barbs, u_barbs, v_barbs, xloc=1.0)
b.set_zorder(4)
skew.shade_cape(prs0, T0-273.15, T0_parcel.magnitude[:]-273.15)
skew.shade_cin(prs0*units.hPa, T0*units.K, T0_parcel, Td0*units.K)

skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 30)
skew.ax.set_xlabel("Temperature ($^{\circ}$C)", fontsize=14)
skew.ax.set_ylabel("Pressure (hPa)", fontsize=14)
skew.ax.tick_params(axis='both', which='major', labelsize=14)
# skew.ax.set_title("ERA5 20250820 T21Z (50.75$^{\circ}$,-114.0$^{\circ}$)", fontsize=16)

# ax_hod = inset_axes(skew.ax, '42%', '42%', loc=1)
# H = Hodograph(ax_hod, component_range=50.)
# H.add_grid(increment=10)
# hod = H.plot_colormapped(u0, v0, z, cmap='HomeyerRainbow', linewidth=1.5)
# H.ax.set_xlim([-20,50])
# H.ax.set_ylim([-20,50])

iz = np.argmin(abs(zh-12))+1

# ax_hod = skew.ax.inset_axes([0.53, 0.56, 0.42, 0.42])
# ax_cb = skew.ax.inset_axes([1.0, 0.56, 0.025, 0.42])
# ax_hod = skew.ax.inset_axes([0.61, 0.65, 0.35, 0.35], zorder=2)
# ax_cb = skew.ax.inset_axes([1.01, 0.675, 0.024, 0.305])
ax_hod = skew.ax.inset_axes([0.52, 0.64, 0.35, 0.35], zorder=2)
H = Hodograph(ax_hod, component_range=50.)
H.add_grid(increment=10)
hod = H.plot_colormapped(u0[:iz], v0[:iz], zh[:iz], cmap='HomeyerRainbow', linewidth=3)

ax_hod.scatter(lm.magnitude[0], lm.magnitude[1], s=12, color='k', zorder=3)
ax_hod.scatter(rm.magnitude[0], rm.magnitude[1], s=12, color='k', zorder=3)
# ax_hod.scatter(mm.magnitude[0], mm.magnitude[1], s=12, color='k', zorder=3)
ax_hod.text(lm.magnitude[0], lm.magnitude[1], 'LM ', ha='right', va='center_baseline', fontsize=10, fontweight='bold')
ax_hod.text(rm.magnitude[0], rm.magnitude[1], ' RM ', ha='left', va='center_baseline', fontsize=10, fontweight='bold')
# ax_hod.text(mm.magnitude[0], mm.magnitude[1], 'MW ', ha='right', va='center_baseline', fontsize=10, fontweight='bold')

H.ax.set_xlim([-20,40]) # -30,30 for domain relative without umove/vmove
H.ax.set_ylim([-20,40])
# H.ax.set_xticks(np.linspace(-30,30,7))
H.ax.xaxis.set_major_locator(MultipleLocator(20))
H.ax.xaxis.set_minor_locator(MultipleLocator(10))
H.ax.yaxis.set_major_locator(MultipleLocator(20))
H.ax.yaxis.set_minor_locator(MultipleLocator(10))
H.ax.set_xlabel("u (m s$^{-1}$)", fontsize=12)
H.ax.set_ylabel("v (m s$^{-1}$)", fontsize=12)
H.ax.tick_params(axis='both', which='major', labelsize=11)

# cb = plt.colorbar(hod, cax=ax_cb)
divider = make_axes_locatable(ax_hod)
cax = divider.append_axes('top', size="5%", pad=0.0)
cb = plt.colorbar(hod, cax=cax, orientation='horizontal', location='top')

cb.set_label("Height AGL (km)", fontsize=12)
cb.set_ticks(np.linspace(0,12,7))
cb.ax.tick_params(labelsize=12)


str1 = '\n'.join((
    "MLCAPE:               %.0f" % (mlcape.magnitude, ),
    "MLCIN:                  %.0f" % (mlcin.magnitude, ),
    "SBCAPE:               %.0f" % (sbcape.magnitude, ),
    "SBCIN:                  %.0f" % (sbcin.magnitude, ),
    "0-1 km shear:       %.1f" % (s01,),
    "0-3 km shear:     %.1f" % (s03,),
    "0-6 km shear:     %.1f" % (s06,),
    "0-1 km SRH:          %.0f" % (srh01,),
    "0-3 km SRH:        %.0f" % (srh03,),
    "LCL height:          %.0f" % (lcl_height*1000,)
    ))
props = dict(boxstyle='square', facecolor='w', edgecolor='k')
skew.ax.text(-39, 950, str1, fontsize=11, bbox=props)




if figsave:
    plt.savefig("D:/brooks/basestate_domain_relative.png", dpi=300)






















