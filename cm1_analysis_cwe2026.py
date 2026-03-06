# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 12:06:45 2026

@author: mschne28
"""

from CM1utils import *

#%% Overview plot for all 3 simulations - dbz and thrpert


fn = np.linspace(5,37,5)

fp1 = 'C:/Users/mschne28/Documents/cm1out/cwe/freeslip_wk_250m/'
fp2 = 'C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_250m/'
fp3 = 'C:/Users/mschne28/Documents/cm1out/cwe/noslip_wk_250m/'

# if 'semislip' in fp:
#     bbc = 'Semi-slip'
#     sim = 'SEMISLIP'
# elif 'freeslip' in fp:
#     bbc = 'Free-slip'
#     sim = 'FREESLIP'
# elif 'noslip' in fp:
#     bbc = 'No-slip'
#     sim = 'NOSLIP'

# titlestr = f"{bbc}, WK profile, dx=250m"
# titlestr = "New P3 -- Fir, modded"


figsave = False

plot_dbz = False
plot_thr = True


if plot_dbz:
    # fig,ax = plt.subplots(3, 4, figsize=(9.5,7), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
    fig,ax = plt.subplots(3, 5, figsize=(11.75,7), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
    
    for f in fn:
        ### Free-slip ###
        ds = nc.Dataset(fp1 + f"cm1out_{f:06.0f}.nc")
        time = ds.variables['time'][:].data[0]
        xh = ds.variables['xh'][:].data
        yh = ds.variables['yh'][:].data
        zh = ds.variables['zh'][:].data
        iz1 = np.where(zh>1)[0][0]
        iz2 = np.where(zh>2)[0][0]
        iz3 = np.where(zh>3)[0][0]
        
        dbz = ds.variables['dbz'][:].data[0,0,:,:]
        winterp = ds.variables['winterp'][:].data[0,0:iz2,:,:]
        umove = ds.variables['umove'][:].data[0]
        vmove = ds.variables['vmove'][:].data[0]
        ds.close()
        
        xl = [-150,150]
        yl = [-150,150]
        
        
        n = int((f-fn[0])/(fn[1]-fn[0]))
        # i = int(np.floor(n/4))
        # j = int(np.mod(n,4))
        xn = xh + (umove*n*7200/1000)
        yn = yh + (vmove*n*7200/1000)
        
        if f == fn[-1]:
            cb_flag = True
        else:
            cb_flag = False
        
        
        
        plot_contourf(xh, yh, np.ma.masked_array(dbz, dbz<0.1), 'dbz', ax[0,n], levels=np.linspace(0,70,15),
                      datalims=[0,70], xlims=xl, ylims=yl, cmap='HomeyerRainbow', cbar=cb_flag)
        ax[0,n].contour(xh, yh, np.max(winterp, axis=0), levels=[5,10], colors=['dimgray','k'], linestyles='-', linewidths=[0.75,0.75])
        if n == 0:
            l1, = ax[0,0].plot([190,200], [190,200], color='dimgray', linewidth=0.75)
            l2, = ax[0,0].plot([190,200], [190,200], '-k', linewidth=0.75)
            ax[0,0].legend(handles=[l1,l2], labels=['w=5 m/s','w=10 m/s'], loc='upper right', fontsize=10)
        ax[0,n].set_title(f"t = {time:.0f} s")
        if n == 0:
            ax[0,n].set_ylabel('y (km)', fontsize=12)
        
        
        ### Semi-slip
        ds = nc.Dataset(fp2 + f"cm1out_{f:06.0f}.nc")
        dbz = ds.variables['dbz'][:].data[0,0,:,:]
        winterp = ds.variables['winterp'][:].data[0,0:iz2,:,:]
        ds.close()
        
        plot_contourf(xh, yh, np.ma.masked_array(dbz, dbz<0.1), 'dbz', ax[1,n], levels=np.linspace(0,70,15),
                      datalims=[0,70], xlims=xl, ylims=yl, cmap='HomeyerRainbow', cbar=cb_flag)
        ax[1,n].contour(xh, yh, np.max(winterp, axis=0), levels=[5,10], colors=['dimgray','k'], linestyles='-', linewidths=[0.75,0.75])
        if n == 0:
            ax[1,n].set_ylabel('y (km)', fontsize=12)
        
        
        ### No-slip
        ds = nc.Dataset(fp3 + f"cm1out_{f:06.0f}.nc")
        dbz = ds.variables['dbz'][:].data[0,0,:,:]
        winterp = ds.variables['winterp'][:].data[0,0:iz2,:,:]
        ds.close()
        
        plot_contourf(xh, yh, np.ma.masked_array(dbz, dbz<0.1), 'dbz', ax[2,n], levels=np.linspace(0,70,15),
                      datalims=[0,70], xlims=xl, ylims=yl, cmap='HomeyerRainbow', cbar=cb_flag)
        ax[2,n].contour(xh, yh, np.max(winterp, axis=0), levels=[5,10], colors=['dimgray','k'], linestyles='-', linewidths=[0.75,0.75])
        ax[2,n].set_xlabel('x (km)', fontsize=12)
        if n == 0:
            ax[2,n].set_ylabel('y (km)', fontsize=12)
        if (f==fn[-1]) & (figsave):
            fig.savefig(fp2+f"figs/dbz_all.png", dpi=300)
        


if plot_thr:
    fig,ax = plt.subplots(3, 5, figsize=(12,7), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
    
    for f in fn:
        ### Free-slip
        ds = nc.Dataset(fp1 + f"cm1out_{f:06.0f}.nc")
        time = ds.variables['time'][:].data[0]
        xh = ds.variables['xh'][:].data
        yh = ds.variables['yh'][:].data
        zh = ds.variables['zh'][:].data
        iz1 = np.where(zh>1)[0][0]
        iz2 = np.where(zh>2)[0][0]
        iz3 = np.where(zh>3)[0][0]
        
        uinterp = ds.variables['uinterp'][:].data[0,0,:,:]
        vinterp = ds.variables['vinterp'][:].data[0,0,:,:]
        zvort = ds.variables['zvort'][:].data[0,0:iz1,:,:]
        u_gr = uinterp + ds.variables['umove'][:].data[0]
        v_gr = vinterp + ds.variables['vmove'][:].data[0]
        umove = ds.variables['umove'][:].data[0]
        vmove = ds.variables['vmove'][:].data[0]
        thrpert = ds.variables['th'][:].data[0,0,:,:] - ds.variables['th0'][:].data[0,0,:,:]
        # # P3 3-moment scheme
        # if 'qi1' in list(ds.variables.keys()):
        #     thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] - 
        #                 (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] + 
        #                  ds.variables['qi1'][:].data[0,0,:,:] +
        #                  ds.variables['qi2'][:].data[0,0,:,:] + 
        #                  ds.variables['qi3'][:].data[0,0,:,:]))
        #                  # ds.variables['qi4'][:].data[0,0,:,:]))
        # # NSSL 3-moment scheme
        # elif 'qg' in list(ds.variables.keys()):
        #     thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] -
        #                 (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] +
        #                  ds.variables['qi'][:].data[0,0,:,:] + ds.variables['qs'][:].data[0,0,:,:] +
        #                  ds.variables['qg'][:].data[0,0,:,:])) # + ds.variables['qhl'][:].data[0,0,:,:]))
        # thr0 = ds.variables['th0'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,:,:])
        # thrpert = thr - thr0
        # del thr,thr0
        
        ds.close()
        
        
        xl = [-150,150]
        yl = [-150,150]
        
        
        n = int((f-fn[0])/(fn[1]-fn[0]))
        # i = int(np.floor(n/4))
        # j = int(np.mod(n,4))
        xn = xh + (umove*n*7200/1000)
        yn = yh + (vmove*n*7200/1000)
        
        if f == fn[-1]:
            cb_flag = True
        else:
            cb_flag = False
        
        
        qix = 60
        
        plot_contourf(xh, yh, thrpert, 'thpert', ax[0,n], levels=np.linspace(-12,12,25),
                      datalims=[-12,12], xlims=xl, ylims=yl, cmap='balance', cbar=cb_flag)
        ax[0,n].contour(xh, yh, np.max(zvort, axis=0), levels=[0.03], colors='r', linestyles='-', linewidths=1)
        ax[0,n].quiver(xh[::qix], yh[::qix], u_gr[::qix,::qix], v_gr[::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
        if n == 0:
            l3, = ax[0,0].plot([190,200], [190,200], '-r', linewidth=1)
            ax[0,0].legend(handles=[l3], labels=["\u03B6=0.03 s$^{-1}$"], loc='upper right', fontsize=10)
        ax[0,n].set_title(f"t = {time:.0f} s")
        if n == 0:
            ax[0,n].set_ylabel('y (km)', fontsize=12)
        
        
        
        
        ### Semi-slip
        ds = nc.Dataset(fp2 + f"cm1out_{f:06.0f}.nc")
        uinterp = ds.variables['uinterp'][:].data[0,0,:,:]
        vinterp = ds.variables['vinterp'][:].data[0,0,:,:]
        u_gr = uinterp + ds.variables['umove'][:].data[0]
        v_gr = vinterp + ds.variables['vmove'][:].data[0]
        zvort = ds.variables['zvort'][:].data[0,0:iz1,:,:]
        thrpert = ds.variables['th'][:].data[0,0,:,:] - ds.variables['th0'][:].data[0,0,:,:]
        # # P3 3-moment scheme
        # if 'qi1' in list(ds.variables.keys()):
        #     thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] - 
        #                 (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] + 
        #                  ds.variables['qi1'][:].data[0,0,:,:] +
        #                  ds.variables['qi2'][:].data[0,0,:,:] + 
        #                  ds.variables['qi3'][:].data[0,0,:,:]))
        #                  # ds.variables['qi4'][:].data[0,0,:,:]))
        # # NSSL 3-moment scheme
        # elif 'qg' in list(ds.variables.keys()):
        #     thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] -
        #                 (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] +
        #                  ds.variables['qi'][:].data[0,0,:,:] + ds.variables['qs'][:].data[0,0,:,:] +
        #                  ds.variables['qg'][:].data[0,0,:,:])) # + ds.variables['qhl'][:].data[0,0,:,:]))
        # thr0 = ds.variables['th0'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,:,:])
        # thrpert = thr - thr0
        # del thr,thr0
        ds.close()
        
        
        plot_contourf(xh, yh, thrpert, 'thpert', ax[1,n], levels=np.linspace(-12,12,25),
                      datalims=[-12,12], xlims=xl, ylims=yl, cmap='balance', cbar=cb_flag)
        ax[1,n].contour(xh, yh, np.max(zvort, axis=0), levels=[0.03], colors='r', linestyles='-', linewidths=1)
        ax[1,n].quiver(xh[::qix], yh[::qix], u_gr[::qix,::qix], v_gr[::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
        if n == 0:
            ax[1,n].set_ylabel('y (km)', fontsize=12)
        
        
        
        ### No-slip
        ds = nc.Dataset(fp3 + f"cm1out_{f:06.0f}.nc")
        uinterp = ds.variables['uinterp'][:].data[0,1,:,:]
        vinterp = ds.variables['vinterp'][:].data[0,1,:,:]
        zvort = ds.variables['zvort'][:].data[0,0:iz1,:,:]
        u_gr = uinterp + ds.variables['umove'][:].data[0]
        v_gr = vinterp + ds.variables['vmove'][:].data[0]
        thrpert = ds.variables['th'][:].data[0,0,:,:] - ds.variables['th0'][:].data[0,0,:,:]
        # # P3 3-moment scheme
        # if 'qi1' in list(ds.variables.keys()):
        #     thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] - 
        #                 (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] + 
        #                  ds.variables['qi1'][:].data[0,0,:,:] +
        #                  ds.variables['qi2'][:].data[0,0,:,:] + 
        #                  ds.variables['qi3'][:].data[0,0,:,:]))
        #                  # ds.variables['qi4'][:].data[0,0,:,:]))
        # # NSSL 3-moment scheme
        # elif 'qg' in list(ds.variables.keys()):
        #     thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] -
        #                 (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] +
        #                  ds.variables['qi'][:].data[0,0,:,:] + ds.variables['qs'][:].data[0,0,:,:] +
        #                  ds.variables['qg'][:].data[0,0,:,:])) # + ds.variables['qhl'][:].data[0,0,:,:]))
        # thr0 = ds.variables['th0'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,:,:])
        # thrpert = thr - thr0
        # del thr,thr0
        ds.close()
        
        
        plot_contourf(xh, yh, thrpert, 'thpert', ax[2,n], levels=np.linspace(-12,12,25),
                      datalims=[-12,12], xlims=xl, ylims=yl, cmap='balance', cbar=cb_flag)
        ax[2,n].contour(xh, yh, np.max(zvort, axis=0), levels=[0.03], colors='r', linestyles='-', linewidths=1)
        ax[2,n].quiver(xh[::qix], yh[::qix], u_gr[::qix,::qix], v_gr[::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
        ax[2,n].set_xlabel('x (km)', fontsize=12)
        if n == 0:
            ax[2,n].set_ylabel('y (km)', fontsize=12)
        if (f==fn[-1]) & (figsave):
            fig.savefig(fp2+f"figs/thrpert_all.png", dpi=300)


#%% Time series plots

fp1 = 'C:/Users/mschne28/Documents/cm1out/cwe/freeslip_wk_250m/'
ds = nc.Dataset(fp1+f"cm1out_stats.nc")
time = ds.variables['mtime'][:].data
wmax500_fs = ds.variables['wmax500'][:].data #max w at 500 m
wmax1000_fs = ds.variables['wmax1000'][:].data #max w at 1000 m
wmax2500_fs = ds.variables['wmax2500'][:].data #max w at 2500 m
wmax5000_fs = ds.variables['wmax5000'][:].data #max w at 5000 m
swspmax_fs = ds.variables['swspmax'][:].data #max sfc wspd
vortsfc_fs = ds.variables['vortsfc'][:].data #max sfc vort
vort1km_fs = ds.variables['vort1km'][:].data #max 1km vort
vort2km_fs = ds.variables['vort2km'][:].data #max 2km vort
vort3km_fs = ds.variables['vort3km'][:].data #max 3km vort
sthpmin_fs = ds.variables['sthpmin'][:].data #min sfc thpert
pratemax_fs = ds.variables['pratemax'][:].data #max sfc rain rate
sratemax_fs = ds.variables['sratemax'][:].data #max sfc hail rate
ds.close()

stdev_fs = {'wmax500':np.std(wmax500_fs), 'wmax1000':np.std(wmax1000_fs), 'wmax2500':np.std(wmax2500_fs), 'wmax5000':np.std(wmax5000_fs),
          'vortsfc':np.std(vortsfc_fs), 'vort1km':np.std(vort1km_fs), 'vort2km':np.std(vort2km_fs), 'vort3km':np.std(vort3km_fs),
          'swspmax':np.std(swspmax_fs), 'sthpmin':np.std(sthpmin_fs), 'pratemax':np.std(pratemax_fs), 'sratemax':np.std(sratemax_fs)}
var_fs = {'wmax500':np.var(wmax500_fs), 'wmax1000':np.var(wmax1000_fs), 'wmax2500':np.var(wmax2500_fs), 'wmax5000':np.var(wmax5000_fs),
          'vortsfc':np.var(vortsfc_fs), 'vort1km':np.var(vort1km_fs), 'vort2km':np.var(vort2km_fs), 'vort3km':np.var(vort3km_fs),
          'swspmax':np.var(swspmax_fs), 'sthpmin':np.var(sthpmin_fs), 'pratemax':np.var(pratemax_fs), 'sratemax':np.var(sratemax_fs)}


fp2 = 'C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_250m/'
ds = nc.Dataset(fp2+f"cm1out_stats.nc")
wmax500_ss = ds.variables['wmax500'][:].data #max w at 500 m
wmax1000_ss = ds.variables['wmax1000'][:].data #max w at 1000 m
wmax2500_ss = ds.variables['wmax2500'][:].data #max w at 2500 m
wmax5000_ss = ds.variables['wmax5000'][:].data #max w at 5000 m
swspmax_ss = ds.variables['swspmax'][:].data #max sfc wspd
vortsfc_ss = ds.variables['vortsfc'][:].data #max sfc vort
vort1km_ss = ds.variables['vort1km'][:].data #max 1km vort
vort2km_ss = ds.variables['vort2km'][:].data #max 2km vort
vort3km_ss = ds.variables['vort3km'][:].data #max 3km vort
sthpmin_ss = ds.variables['sthpmin'][:].data #min sfc thpert
pratemax_ss = ds.variables['pratemax'][:].data #max sfc rain rate
sratemax_ss = ds.variables['sratemax'][:].data #max sfc hail rate
ds.close()

stdev_ss = {'wmax500':np.std(wmax500_ss), 'wmax1000':np.std(wmax1000_ss), 'wmax2500':np.std(wmax2500_ss), 'wmax5000':np.std(wmax5000_ss),
          'vortsfc':np.std(vortsfc_ss), 'vort1km':np.std(vort1km_ss), 'vort2km':np.std(vort2km_ss), 'vort3km':np.std(vort3km_ss),
          'swspmax':np.std(swspmax_ss), 'sthpmin':np.std(sthpmin_ss), 'pratemax':np.std(pratemax_ss), 'sratemax':np.std(sratemax_ss)}
var_ss = {'wmax500':np.var(wmax500_ss), 'wmax1000':np.var(wmax1000_ss), 'wmax2500':np.var(wmax2500_ss), 'wmax5000':np.var(wmax5000_ss),
          'vortsfc':np.var(vortsfc_ss), 'vort1km':np.var(vort1km_ss), 'vort2km':np.var(vort2km_ss), 'vort3km':np.var(vort3km_ss),
          'swspmax':np.var(swspmax_ss), 'sthpmin':np.var(sthpmin_ss), 'pratemax':np.var(pratemax_ss), 'sratemax':np.var(sratemax_ss)}


fp3 = 'C:/Users/mschne28/Documents/cm1out/cwe/noslip_wk_250m/'
ds = nc.Dataset(fp3+f"cm1out_stats.nc")
wmax500_ns = ds.variables['wmax500'][:].data #max w at 500 m
wmax1000_ns = ds.variables['wmax1000'][:].data #max w at 1000 m
wmax2500_ns = ds.variables['wmax2500'][:].data #max w at 2500 m
wmax5000_ns = ds.variables['wmax5000'][:].data #max w at 5000 m
swspmax_ns = ds.variables['swspmax'][:].data #max sfc wspd
vortsfc_ns = ds.variables['vortsfc'][:].data #max sfc vort
vort1km_ns = ds.variables['vort1km'][:].data #max 1km vort
vort2km_ns = ds.variables['vort2km'][:].data #max 2km vort
vort3km_ns = ds.variables['vort3km'][:].data #max 3km vort
sthpmin_ns = ds.variables['sthpmin'][:].data #min sfc thpert
pratemax_ns = ds.variables['pratemax'][:].data #max sfc rain rate
sratemax_ns = ds.variables['sratemax'][:].data #max sfc hail rate
ds.close()

stdev_ns = {'wmax500':np.std(wmax500_ns), 'wmax1000':np.std(wmax1000_ns), 'wmax2500':np.std(wmax2500_ns), 'wmax5000':np.std(wmax5000_ns),
          'vortsfc':np.std(vortsfc_ns), 'vort1km':np.std(vort1km_ns), 'vort2km':np.std(vort2km_ns), 'vort3km':np.std(vort3km_ns),
          'swspmax':np.std(swspmax_ns), 'sthpmin':np.std(sthpmin_ns), 'pratemax':np.std(pratemax_ns), 'sratemax':np.std(sratemax_ns)}
var_ns = {'wmax500':np.var(wmax500_ns), 'wmax1000':np.var(wmax1000_ns), 'wmax2500':np.var(wmax2500_ns), 'wmax5000':np.var(wmax5000_ns),
          'vortsfc':np.var(vortsfc_ns), 'vort1km':np.var(vort1km_ns), 'vort2km':np.var(vort2km_ns), 'vort3km':np.var(vort3km_ns),
          'swspmax':np.var(swspmax_ns), 'sthpmin':np.var(sthpmin_ns), 'pratemax':np.var(pratemax_ns), 'sratemax':np.var(sratemax_ns)}



figsave = False


### Vorticity time series

fig,ax = plt.subplots(3, 1, figsize=(10,9), sharex=True, layout='constrained')

l1,= ax[0].plot(time, movmean(vortsfc_fs,5), 'k', linewidth=2)
l2,= ax[0].plot(time, movmean(vortsfc_ns,5), 'dodgerblue', linewidth=2)
l3,= ax[0].plot(time, movmean(vortsfc_ss,5), 'crimson', linewidth=2)
# l1,= ax[0].plot(time, vortsfc_fs, 'k', linewidth=2)
# l2,= ax[0].plot(time, vortsfc_ns, 'dodgerblue', linewidth=2)
# l3,= ax[0].plot(time, vortsfc_ss, 'crimson', linewidth=2)
ax[0].set_xlim([0,32400])
ax[0].set_ylim([0,0.2])
# ax[0].set_xlabel('Time (s)', fontsize=14)
ax[0].set_ylabel("Vorticity (1/s)", fontsize=14)
ax[0].tick_params(axis='both', labelsize=12)
ax[0].set_title(f"Max. 10-m vertical vorticity", fontsize=16)
ax[0].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[0].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[0].xaxis.set_major_locator(MultipleLocator(3600))
ax[0].xaxis.set_minor_locator(MultipleLocator(900))
ax[0].yaxis.set_major_locator(MultipleLocator(0.05))
ax[0].yaxis.set_minor_locator(MultipleLocator(0.025))
ax[0].legend(handles=[l1,l2,l3], labels=['FREESLIP','NOSLIP','SEMISLIP'],
             loc='upper left', fontsize=14)

l4,= ax[1].plot(time, movmean(vort1km_fs,5), 'k', linewidth=2)
l5,= ax[1].plot(time, movmean(vort1km_ns,5), 'dodgerblue', linewidth=2)
l6,= ax[1].plot(time, movmean(vort1km_ss,5), 'crimson', linewidth=2)
# l4,= ax[1].plot(time, vort1km_fs, 'k', linewidth=2)
# l5,= ax[1].plot(time, vort1km_ns, 'dodgerblue', linewidth=2)
# l6,= ax[1].plot(time, vort1km_ss, 'crimson', linewidth=2)
ax[1].set_xlim([0,32400])
ax[1].set_ylim([0,0.2])
# ax[1].set_xlabel('Time (s)', fontsize=14)
ax[1].set_ylabel("Vorticity (1/s)", fontsize=14)
ax[1].tick_params(axis='both', labelsize=12)
ax[1].set_title(f"Max. 1-km vertical vorticity", fontsize=16)
ax[1].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[1].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[1].xaxis.set_major_locator(MultipleLocator(3600))
ax[1].xaxis.set_minor_locator(MultipleLocator(900))
ax[1].yaxis.set_major_locator(MultipleLocator(0.05))
ax[1].yaxis.set_minor_locator(MultipleLocator(0.025))
# ax[1].legend(handles=[l4,l5,l6], labels=['FREESLIP','NOSLIP','SEMISLIP'],
#              loc='upper left', fontsize=14)

l7,= ax[2].plot(time, movmean(vort3km_fs,5), 'k', linewidth=2)
l8,= ax[2].plot(time, movmean(vort3km_ns,5), 'dodgerblue', linewidth=2)
l9,= ax[2].plot(time, movmean(vort3km_ss,5), 'crimson', linewidth=2)
# l7,= ax[2].plot(time, vort3km_fs, 'k', linewidth=2)
# l8,= ax[2].plot(time, vort3km_ns, 'dodgerblue', linewidth=2)
# l9,= ax[2].plot(time, vort3km_ss, 'crimson', linewidth=2)
ax[2].set_xlim([0,32400])
ax[2].set_ylim([0,0.2])
ax[2].set_xlabel('Time (s)', fontsize=14)
ax[2].set_ylabel("Vorticity (1/s)", fontsize=14)
ax[2].tick_params(axis='both', labelsize=12)
ax[2].set_title(f"Max. 3-km vertical vorticity", fontsize=16)
ax[2].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[2].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[2].xaxis.set_major_locator(MultipleLocator(3600))
ax[2].xaxis.set_minor_locator(MultipleLocator(900))
ax[2].yaxis.set_major_locator(MultipleLocator(0.05))
ax[2].yaxis.set_minor_locator(MultipleLocator(0.025))
# ax[2].legend(handles=[l7,l8,l9], labels=['FREESLIP','NOSLIP','SEMISLIP'],
#              loc='upper left', fontsize=14)

if figsave:
    plt.savefig(fp2+'figs/zeta_all_timeseries.png', dpi=300)
# plt.show()



### Updraft time series

fig,ax = plt.subplots(3, 1, figsize=(10,9), sharex=True, layout='constrained')

l1,= ax[0].plot(time, movmean(wmax1000_fs,5), 'k', linewidth=2)
l2,= ax[0].plot(time, movmean(wmax1000_ns,5), 'dodgerblue', linewidth=2)
l3,= ax[0].plot(time, movmean(wmax1000_ss,5), 'crimson', linewidth=2)
ax[0].set_xlim([0,32400])
ax[0].set_ylim([0,30])
# ax[0].set_xlabel('Time (s)', fontsize=14)
ax[0].set_ylabel("w (m/s)", fontsize=14)
ax[0].tick_params(axis='both', labelsize=12)
ax[0].set_title(f"Max. 1-km updraft speed", fontsize=16)
ax[0].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[0].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[0].xaxis.set_major_locator(MultipleLocator(3600))
ax[0].xaxis.set_minor_locator(MultipleLocator(900))
ax[0].yaxis.set_major_locator(MultipleLocator(5))
ax[0].yaxis.set_minor_locator(MultipleLocator(2.5))
ax[0].legend(handles=[l1,l2,l3], labels=['FREESLIP','NOSLIP','SEMISLIP'],
             loc='upper left', fontsize=14)

l4,= ax[1].plot(time, movmean(wmax2500_fs,5), 'k', linewidth=2)
l5,= ax[1].plot(time, movmean(wmax2500_ns,5), 'dodgerblue', linewidth=2)
l6,= ax[1].plot(time, movmean(wmax2500_ss,5), 'crimson', linewidth=2)
ax[1].set_xlim([0,32400])
ax[1].set_ylim([10,40])
# ax[1].set_xlabel('Time (s)', fontsize=14)
ax[1].set_ylabel("w (m/s)", fontsize=14)
ax[1].tick_params(axis='both', labelsize=12)
ax[1].set_title(f"Max. 2.5-km updraft speed", fontsize=16)
ax[1].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[1].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[1].xaxis.set_major_locator(MultipleLocator(3600))
ax[1].xaxis.set_minor_locator(MultipleLocator(900))
ax[1].yaxis.set_major_locator(MultipleLocator(5))
ax[1].yaxis.set_minor_locator(MultipleLocator(2.5))
# ax[1].legend(handles=[l4,l5,l6], labels=['FREESLIP','NOSLIP','SEMISLIP'],
#              loc='upper left', fontsize=14)

l7,= ax[2].plot(time, movmean(wmax5000_fs,5), 'k', linewidth=2)
l8,= ax[2].plot(time, movmean(wmax5000_ns,5), 'dodgerblue', linewidth=2)
l9,= ax[2].plot(time, movmean(wmax5000_ss,5), 'crimson', linewidth=2)
ax[2].set_xlim([0,32400])
ax[2].set_ylim([20,50]) #[20,50] for 5 km
ax[2].set_xlabel('Time (s)', fontsize=14)
ax[2].set_ylabel("w (m/s)", fontsize=14)
ax[2].tick_params(axis='both', labelsize=12)
ax[2].set_title(f"Max. 5-km updraft speed", fontsize=16)
ax[2].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[2].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[2].xaxis.set_major_locator(MultipleLocator(3600))
ax[2].xaxis.set_minor_locator(MultipleLocator(900))
ax[2].yaxis.set_major_locator(MultipleLocator(5))
ax[2].yaxis.set_minor_locator(MultipleLocator(2.5))
# ax[2].legend(handles=[l7,l8,l9], labels=['FREESLIP','NOSLIP','SEMISLIP'],
#              loc='upper left', fontsize=14)

if figsave:
    plt.savefig(fp2+'figs/w_all_timeseries.png', dpi=300)




### Cold pool strength metrics

fig,ax = plt.subplots(3, 1, figsize=(10,9), sharex=True, layout='constrained')

l1,= ax[0].plot(time, movmean(swspmax_fs,5), 'k', linewidth=2)
l2,= ax[0].plot(time, movmean(swspmax_ns,5), 'dodgerblue', linewidth=2)
l3,= ax[0].plot(time, movmean(swspmax_ss,5), 'crimson', linewidth=2)
ax[0].set_xlim([0,32400])
ax[0].set_ylim([0,60])
# ax[0].set_xlabel('Time (s)', fontsize=14)
ax[0].set_ylabel("Wind speed (m/s)", fontsize=14)
ax[0].tick_params(axis='both', labelsize=12)
ax[0].set_title(f"Max. 10-m wind speed", fontsize=16)
ax[0].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[0].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[0].xaxis.set_major_locator(MultipleLocator(3600))
ax[0].xaxis.set_minor_locator(MultipleLocator(900))
ax[0].yaxis.set_major_locator(MultipleLocator(10))
ax[0].yaxis.set_minor_locator(MultipleLocator(5))
ax[0].legend(handles=[l1,l2,l3], labels=['FREESLIP','NOSLIP','SEMISLIP'],
             loc='upper left', fontsize=14)

l4,= ax[1].plot(time, sthpmin_fs, 'k', linewidth=2)
l5,= ax[1].plot(time, sthpmin_ns, 'dodgerblue', linewidth=2)
l6,= ax[1].plot(time, sthpmin_ss, 'crimson', linewidth=2)
ax[1].set_xlim([0,32400])
ax[1].set_ylim([-15,0])
# ax[1].set_xlabel('Time (s)', fontsize=14)
ax[1].set_ylabel("\u03B8' (K)", fontsize=14)
ax[1].tick_params(axis='both', labelsize=12)
ax[1].set_title(f"Min. 10-m temperature perturbation", fontsize=16)
ax[1].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[1].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[1].xaxis.set_major_locator(MultipleLocator(3600))
ax[1].xaxis.set_minor_locator(MultipleLocator(900))
ax[1].yaxis.set_major_locator(MultipleLocator(5))
ax[1].yaxis.set_minor_locator(MultipleLocator(2.5))
# ax[1].legend(handles=[l4,l5,l6], labels=['FREESLIP','NOSLIP','SEMISLIP'],
#              loc='upper left', fontsize=14)

l7,= ax[2].plot(time, movmean(pratemax_fs+sratemax_fs,5), 'k', linewidth=2)
l8,= ax[2].plot(time, movmean(pratemax_ns+sratemax_ns,5), 'dodgerblue', linewidth=2)
l9,= ax[2].plot(time, movmean(pratemax_ss+sratemax_ss,5), 'crimson', linewidth=2)
# l7,= ax[2].plot(time, pratemax_fs, 'k', linewidth=2)
# l8,= ax[2].plot(time, sratemax_fs, '--k', linewidth=2)
# ax[2].plot(time, pratemax_ns, 'dodgerblue', linewidth=2)
# ax[2].plot(time, sratemax_ns, 'dodgerblue, linewidth=2, linestyle='--')
# ax[2].plot(time, pratemax_ss, 'crimson', linewidth=2)
# ax[2].plot(time, sratemax_ss, 'crimson', linewidth=2, linestyle='--')
ax[2].set_xlim([0,32400])
ax[2].set_ylim([0,0.1]) #[20,50] for 5 km
ax[2].set_xlabel('Time (s)', fontsize=14)
ax[2].set_ylabel("Precip. rate (kg m$^{-2}$ s$^{-1}$)", fontsize=14)
ax[2].tick_params(axis='both', labelsize=12)
ax[2].set_title(f"Max. surface precip. rate", fontsize=16)
ax[2].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[2].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[2].xaxis.set_major_locator(MultipleLocator(3600))
ax[2].xaxis.set_minor_locator(MultipleLocator(900))
ax[2].yaxis.set_major_locator(MultipleLocator(0.02))
ax[2].yaxis.set_minor_locator(MultipleLocator(0.01))
# ax[2].legend(handles=[l7,l8,l9], labels=['FREESLIP','NOSLIP','SEMISLIP'],
#              loc='upper left', fontsize=14)
# ax[2].legend(handles=[l7,l8], labels=['Rain','Ice'],
#              loc='upper left', fontsize=14)

if figsave:
    plt.savefig(fp2+'figs/coldpool_all_timeseries.png', dpi=300)




#%% Swaths?

# Freeslip
fp1 = 'C:/Users/mschne28/Documents/cm1out/cwe/freeslip_wk_250m/'
ds = nc.Dataset(fp1+'cm1out_000029.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
sws_fs = ds.variables['sws'][:].data[0,:,:] #max sfc wind
svs_fs = ds.variables['svs'][:].data[0,:,:] #max sfc vort
sus_fs = ds.variables['sus'][:].data[0,:,:] #max 5km updraft
shs_fs = ds.variables['shs'][:].data[0,:,:] #max integrated UH
ds.close()


# Semislip
fp2 = 'C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_250m/'
ds = nc.Dataset(fp2+'cm1out_000029.nc')
sws_ss = ds.variables['sws'][:].data[0,:,:]
svs_ss = ds.variables['svs'][:].data[0,:,:]
sus_ss = ds.variables['sus'][:].data[0,:,:]
shs_ss = ds.variables['shs'][:].data[0,:,:]
ds.close()

# Noslip
fp3 = 'C:/Users/mschne28/Documents/cm1out/cwe/noslip_wk_250m/'
ds = nc.Dataset(fp3+'cm1out_000029.nc')
sws_ns = ds.variables['sws'][:].data[0,:,:]
svs_ns = ds.variables['svs'][:].data[0,:,:]
sus_ns = ds.variables['sus'][:].data[0,:,:]
shs_ns = ds.variables['shs'][:].data[0,:,:]
ds.close()


xl = [-150,150]
yl = [-150,150]


figsave = False



fig,ax = plt.subplots(3, 3, figsize=(11,11.5), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

# Row 1 - freeslip
plot_contourf(xh, yh, np.ma.masked_array(sus_fs, sus_fs<1), 'w', ax[0,0],
              levels=np.linspace(0,40,21), datalims=[0,40],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max', cbar=False)
# ax[0,0].contour(xh, yh, svs_fs, levels=[0.04], colors='k', linewidths=[0.75])
# l1, = ax[0,0].plot([190,200], [190,200], '-k', linewidth=0.75)
# ax[0,0].legend(handles=[l1], labels=["\u03B6=0.04 s$^{-1}$"], loc='upper right', fontsize=12)
ax[0,0].set_title("5-km updraft", fontsize=16)
ax[0,0].set_ylabel('y (km)', fontsize=14)


plot_contourf(xh, yh, np.ma.masked_array(svs_fs, svs_fs<0.001), 'zvort', ax[0,1],
              levels=np.linspace(0,0.05,21), datalims=[0,0.05],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,0.05,11), extend='max', cbar=False)
ax[0,1].contour(xh, yh, sus_fs, levels=[20,35], colors=['dimgray','k'], linewidths=[0.75,1.5])
l2, = ax[0,1].plot([190,200], [190,200], color='dimgray', linewidth=0.75)
l3, = ax[0,1].plot([190,200], [190,200], '-k', linewidth=1.5)
ax[0,1].legend(handles=[l2,l3], labels=['w=20 m/s','w=35 m/s'], loc='upper right', fontsize=12)
ax[0,1].set_title("Surface vorticity", fontsize=16)


plot_contourf(xh, yh, sws_fs, 'wspd', ax[0,2],
              levels=np.linspace(0,40,21), datalims=[0,40], 
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max', cbar=False)
ax[0,2].contour(xh, yh, sws_fs, levels=[26,33], colors=['navy','k'], linewidths=[0.5,1])
l4, = ax[0,2].plot([190,200], [190,200], color='navy', linewidth=0.5)
l5, = ax[0,2].plot([190,200], [190,200], '-k', linewidth=1)
ax[0,2].legend(handles=[l4,l5], labels=['Severe (26 m/s)','Sig. Severe (33 m/s)'], loc='upper right', fontsize=12)
ax[0,2].set_title("Surface wind speed", fontsize=16)



# Row 2 - semislip
plot_contourf(xh, yh, np.ma.masked_array(sus_ss, sus_ss<1), 'w', ax[1,0],
              levels=np.linspace(0,40,21), datalims=[0,40],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max', cbar=False)
# ax[1,0].contour(xh, yh, svs_ss, levels=[0.04], colors='k', linewidths=[0.75])
ax[1,0].set_ylabel('y (km)', fontsize=14)


plot_contourf(xh, yh, np.ma.masked_array(svs_ss, svs_ss<0.001), 'zvort', ax[1,1],
              levels=np.linspace(0,0.05,21), datalims=[0,0.05],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,0.05,11), extend='max', cbar=False)
ax[1,1].contour(xh, yh, sus_ss, levels=[20,35], colors=['dimgray','k'], linewidths=[0.75,1.5])


plot_contourf(xh, yh, sws_ss, 'wspd', ax[1,2],
              levels=np.linspace(0,40,21), datalims=[0,40], 
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max', cbar=False)
ax[1,2].contour(xh, yh, sws_ss, levels=[26,33], colors=['navy','k'], linewidths=[0.5,1])



# Row 3 - noslip
c1 = plot_contourf(xh, yh, np.ma.masked_array(sus_ns, sus_ns<1), 'w', ax[2,0],
              levels=np.linspace(0,40,21), datalims=[0,40],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max', cbar=False)
# ax[2,0].contour(xh, yh, svs_ns, levels=[0.04], colors='k', linewidths=[0.75])
ax[2,0].set_xlabel('x (km)', fontsize=14)
ax[2,0].set_ylabel('y (km)', fontsize=14)
cb1 = plt.colorbar(c1, ax=ax[2,0], location='bottom', extend='max')
cb1.set_label("w (m s$^{-1}$)", fontsize=14)
cb1.set_ticks(np.linspace(0,40,11))


c2 = plot_contourf(xh, yh, np.ma.masked_array(svs_ns, svs_ns<0.001), 'zvort', ax[2,1],
              levels=np.linspace(0,0.05,21), datalims=[0,0.05],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,0.05,11), extend='max', cbar=False)
ax[2,1].contour(xh, yh, sus_ns, levels=[20,35], colors=['dimgray','k'], linewidths=[0.75,1.5])
ax[2,1].set_xlabel('x (km)', fontsize=14)
cb2 = plt.colorbar(c2, ax=ax[2,1], location='bottom', extend='max')
cb2.set_label("\u03B6 (s$^{-1}$)", fontsize=14)
cb2.set_ticks(np.linspace(0,0.05,11))
cb2.formatter.set_powerlimits((0,0))


c3 = plot_contourf(xh, yh, sws_ns, 'wspd', ax[2,2],
              levels=np.linspace(0,40,21), datalims=[0,40], 
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max', cbar=False)
ax[2,2].contour(xh, yh, sws_ns, levels=[26,33], colors=['navy','k'], linewidths=[0.5,1])
ax[2,2].set_xlabel('x (km)', fontsize=14)
cb3 = plt.colorbar(c3, ax=ax[2,2], location='bottom', extend='max')
cb3.set_label("Wind speed (m s$^{-1}$)", fontsize=14)
cb3.set_ticks(np.linspace(0,40,11))


if figsave:
    plt.savefig(fp2+'swaths_all.png', dpi=300)

plt.show()


#%% Other plots?


















