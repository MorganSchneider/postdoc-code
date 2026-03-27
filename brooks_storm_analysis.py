# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 15:05:24 2026

@author: mschne28
"""

from CM1utils import *

#%% Plot simulations


fp = 'C:/Users/mschne28/Documents/cm1out/brooks/era5_125m_test1/'

ds = nc.Dataset(fp+'cm1out_000013.nc')
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
ax.contour(xh, yh, thpert, levels=[-5,-2], colors='dodgerblue', linewidths=[1.5,1], linestyles='--')
ax.contour(xh, yh, wmax, levels=[5,10], colors=['#516572','k'], linewidths=[1,1.25])
ax.contour(xh, yh, vortmax, levels=[0.025,0.05], colors=['maroon','r'], linewidths=[1,1.5], linestyles='-')
ax.set_xlim([-100,100])
ax.set_ylim([-100,100])
ax.set_title(f"Surface dBZ (t = {time/60:.0f} min)", fontsize=12)
ax.set_xlabel("x (km)", fontsize=12)
ax.set_ylabel("y (km)", fontsize=12)

l7,= ax.plot([-150,-149], [-150,-149], color='dodgerblue', linestyle='--', linewidth=1)
l8,= ax.plot([-150,-149], [-150,-149], color='dodgerblue', linestyle='--', linewidth=1.5)
l1,= ax.plot([-150,-149], [-150,-149], color='#516572', linestyle='-', linewidth=1)
l2,= ax.plot([-150,-149], [-150,-149], color='k', linestyle='-', linewidth=1.25)
# l3,= ax.plot([-150,-149], [-150,-149], color='k', linestyle='-', linewidth=2.5)
l4,= ax.plot([-150,-149], [-150,-149], color='maroon', linestyle='-', linewidth=1)
l5,= ax.plot([-150,-149], [-150,-149], color='r', linestyle='-', linewidth=1.5)
# l6,= ax.plot([-150,-149], [-150,-149], color='maroon', linestyle='-', linewidth=1.5)

ax.legend(handles=[l1,l2,l4,l5,l7,l8],
          labels=["w=5 m s$^{-1}$", "w=10 m s$^{-1}$",
                  "\u03B6=0.025 s$^{-1}$", "\u03B6=0.05 s$^{-1}$",
                  "\u03B8'=-2 K", "\u03B8'=-5 K"],
          ncols=3, fontsize=10, loc='lower right')
# plt.savefig(fp+f"dbz_{time/60:.0f}min_125m_run.png", dpi=300)



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
ax.legend(handles=[l1,l2,l3], labels=["\u03B6=0.025 s$^{-1}$", "\u03B6=0.05 s$^{-1}$"],
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


fig,ax = plt.subplots(1, 1, figsize=(10,5), layout='constrained')
l1, = ax.plot(time, vort1km, '-r', linewidth=2)
l2, = ax.plot(time, vort3km, '-b', linewidth=2)
l3, = ax.plot(time, vort5km, '-k', linewidth=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax.xaxis.set_major_locator(MultipleLocator(1800))
ax.xaxis.set_minor_locator(MultipleLocator(300))
ax.yaxis.set_major_locator(MultipleLocator(0.04))
ax.yaxis.set_minor_locator(MultipleLocator(0.01))
ax.set_xlim([0,14400])
ax.set_ylim([0,0.2])
ax.legend(handles=[l1,l2,l3], labels=["\u03B6(1 km)", "\u03B6(3 km)", "\u03B6(5 km)"], loc='lower right', fontsize=12)
ax.set_title("Domain-maximum \u03B6", fontsize=14)
ax.set_xlabel('Time (s)', fontsize=12)
ax.set_ylabel("\u03B6 (1/s)", fontsize=12)
plt.savefig(fp+'vortmax_timeseries.png', dpi=300)


fig,ax = plt.subplots(1, 1, figsize=(10,5), layout='constrained')
l1, = ax.plot(time, wmax1000, '-r', linewidth=2)
l2, = ax.plot(time, wmax2500, '-b', linewidth=2)
l3, = ax.plot(time, wmax5000, '-k', linewidth=2)
ax.grid(visible=True, which='major', color='dimgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax.xaxis.set_major_locator(MultipleLocator(1800))
ax.xaxis.set_minor_locator(MultipleLocator(300))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.set_xlim([0,14400])
ax.set_ylim([0,30])
ax.legend(handles=[l1,l2,l3], labels=['w(1 km)', 'w(2.5 km)', 'w(5 km)'], loc='lower right', fontsize=12)
ax.set_title('Domain-maximum w', fontsize=14)
ax.set_xlabel('Time (s)', fontsize=12)
ax.set_ylabel('w (m/s)', fontsize=12)
plt.savefig(fp+'wmax_timeseries.png', dpi=300)

plt.show()



#%% Check swaths

fp = 'C:/Users/mschne28/Documents/cm1out/brooks/era5_125m_test1/'

ds = nc.Dataset(fp+'cm1out_000017.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
sws = ds.variables['sws'][:].data[0,:,:]
svs = ds.variables['svs'][:].data[0,:,:]
sus = ds.variables['sus'][:].data[0,:,:]
shs = ds.variables['shs2'][:].data[0,:,:]
ds.close()

xl = [-100,100]
yl = [-100,100]


# fig,ax = plt.subplots(1, 1, figsize=(8,6), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))
# plot_contourf(xh, yh, sus, 'w', ax, levels=np.linspace(0,40,21), datalims=[0,40],
#               xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max')
# ax.set_title("5-km updraft swath", fontsize=16)
# ax.set_xlabel('x (km)', fontsize=14)
# ax.set_ylabel('y (km)', fontsize=14)


fig,ax = plt.subplots(1, 1, figsize=(8,6), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))
plot_contourf(xh, yh, shs, 'uh', ax, levels=np.linspace(0,1000,21), datalims=[0,1000],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,1000,11), extend='max')
ax.set_title("Integrated updraft helicity swath (0-4 h)", fontsize=16)
ax.set_xlabel('x (km)', fontsize=14)
ax.set_ylabel('y (km)', fontsize=14)

plt.show()

#%% Translated swaths

fp = 'C:/Users/mschne28/Documents/cm1out/brooks/era5_125m_test1/'

fn = np.linspace(5,17,4)


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




#%% Check updraft nudging

xh = np.linspace(-60000,60000,121)
yh = np.linspace(-60000,60000,121)

xc1 = 0
xc2 = 0
xc3 = 0
xc4 = 0

yc1 = 0
yc2 = -20000
yc3 = 20000
yc4 = 40000

xr = 10000
yr = 10000

alpha = 0.5
wmax = 15
t1 = 900
t2 = 1200

rxr = 1/xr
ryr = 1/yr

w = np.zeros(shape=(20,len(yh),len(xh)), dtype=float)
wa = np.zeros(shape=(len(yh),len(xh)), dtype=float)
w3d = np.zeros(shape=(len(yh),len(xh)), dtype=float)
wten = np.zeros(shape=(20,len(yh),len(xh)), dtype=float)
wndgten = np.zeros(shape=(len(yh),len(xh)), dtype=float)

time = np.linspace(0,t2,t2+1)




dt = 1.0


for t in time:
    gamm = 1.0
    
    if t >= t1:
        gamm = 1 - (t-t1)/(t2-t1)
    
    tem = alpha*gamm
    
    for j in range(len(yh)):
        for i in range(len(xh)):
            beta1 = np.sqrt( ((xh[i]-xc1)*rxr)**2 + ((yh[j]-yc1)*ryr)**2 )
            beta2 = np.sqrt( ((xh[i]-xc2)*rxr)**2 + ((yh[j]-yc2)*ryr)**2 )
            beta3 = np.sqrt( ((xh[i]-xc3)*rxr)**2 + ((yh[j]-yc3)*ryr)**2 )
            beta4 = np.sqrt( ((xh[i]-xc4)*rxr)**2 + ((yh[j]-yc4)*ryr)**2 )
            beta = min([beta1,beta2,beta3,beta4])
            if beta < 1:
                wmag = wmax * (np.cos(0.5*np.pi*beta)**2)
                wndgten[j,i] = tem * max([wmag-wa[j,i], 0])
            else:
                wndgten[j,i] = 0
            
            if t+dt <= t2:
                w3d[j,i] = wa[j,i] + dt*wndgten[j,i]
    
    if t+dt <= t2:
        if np.mod(t,60) == 0:
            n = int(t/60)
            wten[n,:,:] = wndgten[:,:]
            w[n,:,:] = w3d[:,:]
            
    wa[:,:] = w3d[:,:]


#%% I don't get it this is stupid 

fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, w[4,:,:], 'w', ax, levels=np.linspace(0,15,16), datalims=[0,15], cmap='Reds')
ax.set_title("5 min")

fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, w[9,:,:], 'w', ax, levels=np.linspace(0,15,16), datalims=[0,15], cmap='Reds')
ax.set_title("10 min")

fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, w[14,:,:], 'w', ax, levels=np.linspace(0,15,16), datalims=[0,15], cmap='Reds')
ax.set_title("15 min")

fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, w[19,:,:], 'w', ax, levels=np.linspace(0,15,16), datalims=[0,15], cmap='Reds')
ax.set_title("20 min")





fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, wten[0,:,:], 'w', ax, levels=np.linspace(0,7.5,16), datalims=[0,7.5], cmap='Reds')
ax.set_title("1 min")

fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, wten[4,:,:], 'w', ax, levels=np.linspace(0,1e-15,11), datalims=[0,1e-15], cmap='Reds')
ax.set_title("5 min")

fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, wten[14,:,:], 'w', ax, levels=np.linspace(0,1e-15,11), datalims=[0,1e-15], cmap='Reds')
ax.set_title("15 min")

fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, wten[19,:,:], 'w', ax, levels=np.linspace(0,1e-15,11), datalims=[0,1e-15], cmap='Reds')
ax.set_title("20 min")


plt.show()






