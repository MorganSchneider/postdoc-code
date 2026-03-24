# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 12:06:45 2026

@author: mschne28
"""

from CM1utils import *
import metpy.calc as mc
from metpy.units import units

#%% Overview plot for all 3 simulations - dbz, thrpert, laplacian of thrpert


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
plot_thr = False
plot_del2 = True


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
            fig.savefig(fp2+f"figs/dbz_all_v2.png", dpi=300)
        


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
        #                  ds.variables['qg'][:].data[0,0,:,:] + ds.variables['qhl'][:].data[0,0,:,:]))
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
            fig.savefig(fp2+f"figs/thrpert_all_v2.png", dpi=300)



if plot_del2:
    fig,ax = plt.subplots(3, 5, figsize=(12.5,7), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
    
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
        #                  ds.variables['qg'][:].data[0,0,:,:] + ds.variables['qhl'][:].data[0,0,:,:]))
        # thr0 = ds.variables['th0'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,:,:])
        # thrpert = thr - thr0
        # del thr,thr0
        del2 = mc.laplacian(thrpert*units.K, deltas=(250*units.m, 250*units.m))
        del2thp = del2.magnitude
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
        
        c1 = ax[0,n].contourf(xh, yh, del2thp, levels=np.linspace(-5e-5,5e-5,21), vmin=-5e-5, vmax=5e-5, cmap='balance', antialiased=True)
        c1.set_edgecolor('face')
        if cb_flag:
            cb = plt.colorbar(c1, ax=ax[0,n], extend='both')
            cb.set_label("\u25BD$^2$\u03B8' (K m$^{-2}$)", fontsize=11)
            cb.formatter.set_powerlimits((0,0))
            cb.set_ticks(np.linspace(-5e-5,5e-5,11))
        ax[0,n].set_xlim(xl)
        ax[0,n].set_ylim(yl)
        ax[0,n].quiver(xh[::qix], yh[::qix], u_gr[::qix,::qix], v_gr[::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
        ax[0,n].set_title(f"t = {time:.0f} s")
        if n == 0:
            ax[0,n].set_ylabel('y (km)', fontsize=12)
        
        
        
        
        ### Semi-slip
        ds = nc.Dataset(fp2 + f"cm1out_{f:06.0f}.nc")
        uinterp = ds.variables['uinterp'][:].data[0,0,:,:]
        vinterp = ds.variables['vinterp'][:].data[0,0,:,:]
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
        del2 = mc.laplacian(thrpert*units.K, deltas=(250*units.m, 250*units.m))
        del2thp = del2.magnitude
        ds.close()
        
        
        c2 = ax[1,n].contourf(xh, yh, del2thp, levels=np.linspace(-5e-5,5e-5,21), vmin=-5e-5, vmax=5e-5, cmap='balance', antialiased=True)
        c2.set_edgecolor('face')
        if cb_flag:
            cb = plt.colorbar(c2, ax=ax[1,n], extend='both')
            cb.set_label("\u25BD$^2$\u03B8' (K m$^{-2}$)", fontsize=11)
            cb.formatter.set_powerlimits((0,0))
            cb.set_ticks(np.linspace(-5e-5,5e-5,11))
        ax[1,n].set_xlim(xl)
        ax[1,n].set_ylim(yl)
        ax[1,n].quiver(xh[::qix], yh[::qix], u_gr[::qix,::qix], v_gr[::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
        if n == 0:
            ax[1,n].set_ylabel('y (km)', fontsize=12)
        
        
        
        
        ### No-slip
        ds = nc.Dataset(fp3 + f"cm1out_{f:06.0f}.nc")
        uinterp = ds.variables['uinterp'][:].data[0,1,:,:]
        vinterp = ds.variables['vinterp'][:].data[0,1,:,:]
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
        del2 = mc.laplacian(thrpert*units.K, deltas=(250*units.m, 250*units.m))
        del2thp = del2.magnitude
        ds.close()
        
        
        c3 = ax[2,n].contourf(xh, yh, del2thp, levels=np.linspace(-5e-5,5e-5,21), vmin=-5e-5, vmax=5e-5, cmap='balance', antialiased=True)
        c3.set_edgecolor('face')
        if cb_flag:
            cb = plt.colorbar(c3, ax=ax[2,n], extend='both')
            cb.set_label("\u25BD$^2$\u03B8' (K m$^{-2}$)", fontsize=11)
            cb.formatter.set_powerlimits((0,0))
            cb.set_ticks(np.linspace(-5e-5,5e-5,11))
        ax[2,n].set_xlim(xl)
        ax[2,n].set_ylim(yl)
        ax[2,n].quiver(xh[::qix], yh[::qix], u_gr[::qix,::qix], v_gr[::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
        ax[2,n].set_xlabel('x (km)', fontsize=12)
        if n == 0:
            ax[2,n].set_ylabel('y (km)', fontsize=12)
        if (f==fn[-1]) & (figsave):
            fig.savefig(fp2+f"figs/del2thp_all_v2.png", dpi=300)



#%% Time series and statistics

from scipy import stats


fp1 = 'C:/Users/mschne28/Documents/cm1out/cwe/freeslip_wk_250m/'
ds = nc.Dataset(fp1+f"cm1out_stats.nc")
time = ds.variables['mtime'][:].data
wmax500_fs = ds.variables['wmax500'][:].data #max w at 500 m
wmax1000_fs = ds.variables['wmax1000'][:].data #max w at 1000 m
wmax2500_fs = ds.variables['wmax2500'][:].data #max w at 2500 m
wmax5000_fs = ds.variables['wmax5000'][:].data #max w at 5000 m
wmin500_fs = ds.variables['wmin500'][:].data #min w at 500 m
wmin1000_fs = ds.variables['wmin1000'][:].data #min w at 1000 m
wmin2500_fs = ds.variables['wmin2500'][:].data #min w at 2500 m
wmin5000_fs = ds.variables['wmin5000'][:].data #min w at 5000 m
swspmax_fs = ds.variables['swspmax'][:].data #max sfc wspd
vortsfc_fs = ds.variables['vortsfc'][:].data #max sfc vort
vort1km_fs = ds.variables['vort1km'][:].data #max 1km vort
vort2km_fs = ds.variables['vort2km'][:].data #max 2km vort
vort3km_fs = ds.variables['vort3km'][:].data #max 3km vort
sthpmin_fs = ds.variables['sthpmin'][:].data #min sfc thpert
pratemax_fs = ds.variables['pratemax'][:].data #max sfc rain rate
sratemax_fs = ds.variables['sratemax'][:].data #max sfc hail rate
ds.close()

data_fs = {'wmax500':wmax500_fs, 'wmax1000':wmax1000_fs, 'wmax2500':wmax2500_fs, 'wmax5000':wmax5000_fs,
           'wmin500':wmin500_fs, 'wmin1000':wmin1000_fs, 'wmin2500':wmin2500_fs, 'wmin5000':wmin5000_fs,
          'vortsfc':vortsfc_fs, 'vort1km':vort1km_fs, 'vort2km':vort2km_fs, 'vort3km':vort3km_fs,
          'swspmax':swspmax_fs, 'sthpmin':sthpmin_fs, 'pratemax':pratemax_fs, 'sratemax':sratemax_fs}
stdev_fs = {'wmax500':np.std(wmax500_fs), 'wmax1000':np.std(wmax1000_fs), 'wmax2500':np.std(wmax2500_fs), 'wmax5000':np.std(wmax5000_fs),
            'wmin500':np.std(wmin500_fs), 'wmin1000':np.std(wmin1000_fs), 'wmin2500':np.std(wmin2500_fs), 'wmin5000':np.std(wmin5000_fs),
          'vortsfc':np.std(vortsfc_fs), 'vort1km':np.std(vort1km_fs), 'vort2km':np.std(vort2km_fs), 'vort3km':np.std(vort3km_fs),
          'swspmax':np.std(swspmax_fs), 'sthpmin':np.std(sthpmin_fs), 'pratemax':np.std(pratemax_fs), 'sratemax':np.std(sratemax_fs)}
var_fs = {'wmax500':np.var(wmax500_fs), 'wmax1000':np.var(wmax1000_fs), 'wmax2500':np.var(wmax2500_fs), 'wmax5000':np.var(wmax5000_fs),
          'wmin500':np.var(wmin500_fs), 'wmin1000':np.var(wmin1000_fs), 'wmin2500':np.var(wmin2500_fs), 'wmin5000':np.var(wmin5000_fs),
          'vortsfc':np.var(vortsfc_fs), 'vort1km':np.var(vort1km_fs), 'vort2km':np.var(vort2km_fs), 'vort3km':np.var(vort3km_fs),
          'swspmax':np.var(swspmax_fs), 'sthpmin':np.var(sthpmin_fs), 'pratemax':np.var(pratemax_fs), 'sratemax':np.var(sratemax_fs)}


fp2 = 'C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_250m/'
ds = nc.Dataset(fp2+f"cm1out_stats.nc")
wmax500_ss = ds.variables['wmax500'][:].data #max w at 500 m
wmax1000_ss = ds.variables['wmax1000'][:].data #max w at 1000 m
wmax2500_ss = ds.variables['wmax2500'][:].data #max w at 2500 m
wmax5000_ss = ds.variables['wmax5000'][:].data #max w at 5000 m
wmin500_ss = ds.variables['wmin500'][:].data #min w at 500 m
wmin1000_ss = ds.variables['wmin1000'][:].data #min w at 1000 m
wmin2500_ss = ds.variables['wmin2500'][:].data #min w at 2500 m
wmin5000_ss = ds.variables['wmin5000'][:].data #min w at 5000 m
swspmax_ss = ds.variables['swspmax'][:].data #max sfc wspd
vortsfc_ss = ds.variables['vortsfc'][:].data #max sfc vort
vort1km_ss = ds.variables['vort1km'][:].data #max 1km vort
vort2km_ss = ds.variables['vort2km'][:].data #max 2km vort
vort3km_ss = ds.variables['vort3km'][:].data #max 3km vort
sthpmin_ss = ds.variables['sthpmin'][:].data #min sfc thpert
pratemax_ss = ds.variables['pratemax'][:].data #max sfc rain rate
sratemax_ss = ds.variables['sratemax'][:].data #max sfc hail rate
ds.close()

data_ss = {'wmax500':wmax500_ss, 'wmax1000':wmax1000_ss, 'wmax2500':wmax2500_ss, 'wmax5000':wmax5000_ss,
           'wmin500':wmin500_ss, 'wmin1000':wmin1000_ss, 'wmin2500':wmin2500_ss, 'wmin5000':wmin5000_ss,
          'vortsfc':vortsfc_ss, 'vort1km':vort1km_ss, 'vort2km':vort2km_ss, 'vort3km':vort3km_ss,
          'swspmax':swspmax_ss, 'sthpmin':sthpmin_ss, 'pratemax':pratemax_ss, 'sratemax':sratemax_ss}
stdev_ss = {'wmax500':np.std(wmax500_ss), 'wmax1000':np.std(wmax1000_ss), 'wmax2500':np.std(wmax2500_ss), 'wmax5000':np.std(wmax5000_ss),
            'wmin500':np.std(wmin500_ss), 'wmin1000':np.std(wmin1000_ss), 'wmin2500':np.std(wmin2500_ss), 'wmin5000':np.std(wmin5000_ss),
          'vortsfc':np.std(vortsfc_ss), 'vort1km':np.std(vort1km_ss), 'vort2km':np.std(vort2km_ss), 'vort3km':np.std(vort3km_ss),
          'swspmax':np.std(swspmax_ss), 'sthpmin':np.std(sthpmin_ss), 'pratemax':np.std(pratemax_ss), 'sratemax':np.std(sratemax_ss)}
var_ss = {'wmax500':np.var(wmax500_ss), 'wmax1000':np.var(wmax1000_ss), 'wmax2500':np.var(wmax2500_ss), 'wmax5000':np.var(wmax5000_ss),
          'wmin500':np.var(wmin500_ss), 'wmin1000':np.var(wmin1000_ss), 'wmin2500':np.var(wmin2500_ss), 'wmin5000':np.var(wmin5000_ss),
          'vortsfc':np.var(vortsfc_ss), 'vort1km':np.var(vort1km_ss), 'vort2km':np.var(vort2km_ss), 'vort3km':np.var(vort3km_ss),
          'swspmax':np.var(swspmax_ss), 'sthpmin':np.var(sthpmin_ss), 'pratemax':np.var(pratemax_ss), 'sratemax':np.var(sratemax_ss)}


fp3 = 'C:/Users/mschne28/Documents/cm1out/cwe/noslip_wk_250m/'
ds = nc.Dataset(fp3+f"cm1out_stats.nc")
wmax500_ns = ds.variables['wmax500'][:].data #max w at 500 m
wmax1000_ns = ds.variables['wmax1000'][:].data #max w at 1000 m
wmax2500_ns = ds.variables['wmax2500'][:].data #max w at 2500 m
wmax5000_ns = ds.variables['wmax5000'][:].data #max w at 5000 m
wmin500_ns = ds.variables['wmin500'][:].data #min w at 500 m
wmin1000_ns = ds.variables['wmin1000'][:].data #min w at 1000 m
wmin2500_ns = ds.variables['wmin2500'][:].data #min w at 2500 m
wmin5000_ns = ds.variables['wmin5000'][:].data #min w at 5000 m
swspmax_ns = ds.variables['swspmax'][:].data #max sfc wspd
vortsfc_ns = ds.variables['vortsfc'][:].data #max sfc vort
vort1km_ns = ds.variables['vort1km'][:].data #max 1km vort
vort2km_ns = ds.variables['vort2km'][:].data #max 2km vort
vort3km_ns = ds.variables['vort3km'][:].data #max 3km vort
sthpmin_ns = ds.variables['sthpmin'][:].data #min sfc thpert
pratemax_ns = ds.variables['pratemax'][:].data #max sfc rain rate
sratemax_ns = ds.variables['sratemax'][:].data #max sfc hail rate
ds.close()

data_ns = {'wmax500':wmax500_ns, 'wmax1000':wmax1000_ns, 'wmax2500':wmax2500_ns, 'wmax5000':wmax5000_ns,
           'wmin500':wmin500_ns, 'wmin1000':wmin1000_ns, 'wmin2500':wmin2500_ns, 'wmin5000':wmin5000_ns,
          'vortsfc':vortsfc_ns, 'vort1km':vort1km_ns, 'vort2km':vort2km_ns, 'vort3km':vort3km_ns,
          'swspmax':swspmax_ns, 'sthpmin':sthpmin_ns, 'pratemax':pratemax_ns, 'sratemax':sratemax_ns}
stdev_ns = {'wmax500':np.std(wmax500_ns), 'wmax1000':np.std(wmax1000_ns), 'wmax2500':np.std(wmax2500_ns), 'wmax5000':np.std(wmax5000_ns),
            'wmin500':np.std(wmin500_ns), 'wmin1000':np.std(wmin1000_ns), 'wmin2500':np.std(wmin2500_ns), 'wmin5000':np.std(wmin5000_ns),
          'vortsfc':np.std(vortsfc_ns), 'vort1km':np.std(vort1km_ns), 'vort2km':np.std(vort2km_ns), 'vort3km':np.std(vort3km_ns),
          'swspmax':np.std(swspmax_ns), 'sthpmin':np.std(sthpmin_ns), 'pratemax':np.std(pratemax_ns), 'sratemax':np.std(sratemax_ns)}
var_ns = {'wmax500':np.var(wmax500_ns), 'wmax1000':np.var(wmax1000_ns), 'wmax2500':np.var(wmax2500_ns), 'wmax5000':np.var(wmax5000_ns),
          'wmin500':np.var(wmin500_ns), 'wmin1000':np.var(wmin1000_ns), 'wmin2500':np.var(wmin2500_ns), 'wmin5000':np.var(wmin5000_ns),
          'vortsfc':np.var(vortsfc_ns), 'vort1km':np.var(vort1km_ns), 'vort2km':np.var(vort2km_ns), 'vort3km':np.var(vort3km_ns),
          'swspmax':np.var(swspmax_ns), 'sthpmin':np.var(sthpmin_ns), 'pratemax':np.var(pratemax_ns), 'sratemax':np.var(sratemax_ns)}



t_stats = {'ssfs':{}, 'nsfs':{}, 'ssns':{}}
p_vals_ttest = {'ssfs':{}, 'nsfs':{}, 'ssns':{}}
f_stats_anova = {}
p_vals_anova = {}

for key in list(data_fs.keys()):
    # ANOVA test
    f_stat, p_val = stats.f_oneway(data_fs[key], data_ss[key], data_ns[key])
    f_stats_anova.update({key:f_stat})
    p_vals_anova.update({key:p_val})
    
    # Paired t-test between semislip and freeslip
    t_stat, p_val = stats.ttest_ind(data_ss[key], data_fs[key], equal_var=False)
    t_stats['ssfs'].update({key:t_stat})
    p_vals_ttest['ssfs'].update({key:p_val})
    
    # Paired t-test between noslip and freeslip
    t_stat, p_val = stats.ttest_ind(data_ns[key], data_fs[key], equal_var=False)
    t_stats['nsfs'].update({key:t_stat})
    p_vals_ttest['nsfs'].update({key:p_val})
    
    # Paired t-test between semislip and noslip
    t_stat, p_val = stats.ttest_ind(data_ss[key], data_ns[key], equal_var=False)
    t_stats['ssns'].update({key:t_stat})
    p_vals_ttest['ssns'].update({key:p_val})


t_ssfs_hourly = {}
t_nsfs_hourly = {}
t_ssns_hourly = {}
f_anova_hourly = {}
p_ssfs_hourly = {}
p_nsfs_hourly = {}
p_ssns_hourly = {}
p_anova_hourly = {}
isf = {}
inf = {}
isn = {}
ian = {}

for key in list(data_fs.keys()):
    t_ssfs_hourly.update({key:np.zeros(shape=(9,), dtype=float)})
    t_nsfs_hourly.update({key:np.zeros(shape=(9,), dtype=float)})
    t_ssns_hourly.update({key:np.zeros(shape=(9,), dtype=float)})
    f_anova_hourly.update({key:np.zeros(shape=(9,), dtype=float)})
    p_ssfs_hourly.update({key:np.zeros(shape=(9,), dtype=float)})
    p_nsfs_hourly.update({key:np.zeros(shape=(9,), dtype=float)})
    p_ssns_hourly.update({key:np.zeros(shape=(9,), dtype=float)})
    p_anova_hourly.update({key:np.zeros(shape=(9,), dtype=float)})
    isf.update({key:np.zeros(shape=(9,), dtype=float)})
    inf.update({key:np.zeros(shape=(9,), dtype=float)})
    isn.update({key:np.zeros(shape=(9,), dtype=float)})
    ian.update({key:np.zeros(shape=(9,), dtype=float)})
    
    for i in range(9):
        i1 = i*60
        i2 = (i+1)*60 + 1
        
        # Paired t-test between semislip and freeslip
        t_stat, p_val = stats.ttest_ind(data_ss[key][i1:i2], data_fs[key][i1:i2], equal_var=False)
        t_ssfs_hourly[key][i] = t_stat
        p_ssfs_hourly[key][i] = p_val
        if p_val < 0.05:
            isf[key][i] = 1
        
        # Paired t-test between noslip and freeslip
        t_stat, p_val = stats.ttest_ind(data_ns[key][i1:i2], data_fs[key][i1:i2], equal_var=False)
        t_nsfs_hourly[key][i] = t_stat
        p_nsfs_hourly[key][i] = p_val
        if p_val < 0.05:
            inf[key][i] = 1
        
        # Paired t-test between semislip and noslip
        t_stat, p_val = stats.ttest_ind(data_ss[key][i1:i2], data_ns[key][i1:i2], equal_var=False)
        t_ssns_hourly[key][i] = t_stat
        p_ssns_hourly[key][i] = p_val
        if p_val < 0.05:
            isn[key][i] = 1
        
        # ANOVA test
        f_stat, p_val = stats.f_oneway(data_fs[key][i1:i2], data_ss[key][i1:i2], data_ns[key][i1:i2])
        f_anova_hourly[key][i] = f_stat
        p_anova_hourly[key][i] = p_val
        if p_val < 0.05:
            ian[key][i] = 1




#%% Plot time series

figsave = False


### Vorticity time series

fig,ax = plt.subplots(3, 1, figsize=(10,9), sharex=True, layout='constrained')

l1,= ax[0].plot(time[:-2], movmean(vortsfc_fs,5)[:-2], 'k', linewidth=2)
l2,= ax[0].plot(time[:-2], movmean(vortsfc_ns,5)[:-2], 'dodgerblue', linewidth=2)
l3,= ax[0].plot(time[:-2], movmean(vortsfc_ss,5)[:-2], 'crimson', linewidth=2)
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

l4,= ax[1].plot(time[:-2], movmean(vort1km_fs,5)[:-2], 'k', linewidth=2)
l5,= ax[1].plot(time[:-2], movmean(vort1km_ns,5)[:-2], 'dodgerblue', linewidth=2)
l6,= ax[1].plot(time[:-2], movmean(vort1km_ss,5)[:-2], 'crimson', linewidth=2)
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

l7,= ax[2].plot(time[:-2], movmean(vort3km_fs,5)[:-2], 'k', linewidth=2)
l8,= ax[2].plot(time[:-2], movmean(vort3km_ns,5)[:-2], 'dodgerblue', linewidth=2)
l9,= ax[2].plot(time[:-2], movmean(vort3km_ss,5)[:-2], 'crimson', linewidth=2)
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

l1,= ax[0].plot(time[:-2], movmean(wmax1000_fs,5)[:-2], 'k', linewidth=2)
l2,= ax[0].plot(time[:-2], movmean(wmax1000_ns,5)[:-2], 'dodgerblue', linewidth=2)
l3,= ax[0].plot(time[:-2], movmean(wmax1000_ss,5)[:-2], 'crimson', linewidth=2)
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

l4,= ax[1].plot(time[:-2], movmean(wmax2500_fs,5)[:-2], 'k', linewidth=2)
l5,= ax[1].plot(time[:-2], movmean(wmax2500_ns,5)[:-2], 'dodgerblue', linewidth=2)
l6,= ax[1].plot(time[:-2], movmean(wmax2500_ss,5)[:-2], 'crimson', linewidth=2)
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

l7,= ax[2].plot(time[:-2], movmean(wmax5000_fs,5)[:-2], 'k', linewidth=2)
l8,= ax[2].plot(time[:-2], movmean(wmax5000_ns,5)[:-2], 'dodgerblue', linewidth=2)
l9,= ax[2].plot(time[:-2], movmean(wmax5000_ss,5)[:-2], 'crimson', linewidth=2)
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

l1,= ax[0].plot(time[:-2], movmean(swspmax_fs,5)[:-2], 'k', linewidth=2)
l2,= ax[0].plot(time[:-2], movmean(swspmax_ns,5)[:-2], 'dodgerblue', linewidth=2)
l3,= ax[0].plot(time[:-2], movmean(swspmax_ss,5)[:-2], 'crimson', linewidth=2)
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

l4,= ax[1].plot(time[:-2], sthpmin_fs[:-2], 'k', linewidth=2)
l5,= ax[1].plot(time[:-2], sthpmin_ns[:-2], 'dodgerblue', linewidth=2)
l6,= ax[1].plot(time[:-2], sthpmin_ss[:-2], 'crimson', linewidth=2)
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

l7,= ax[2].plot(time[:-2], movmean(pratemax_fs,5)[:-2], 'k', linewidth=2)
ax[2].plot(time[:-2], movmean(pratemax_ns,5)[:-2], 'dodgerblue', linewidth=2)
ax[2].plot(time[:-2], movmean(pratemax_ss,5)[:-2], 'crimson', linewidth=2)
# l7,= ax[2].plot(time, pratemax_fs, 'k', linewidth=2)
# l8,= ax[2].plot(time, sratemax_fs, '--k', linewidth=2)
# ax[2].plot(time, pratemax_ns, 'dodgerblue', linewidth=2)
# ax[2].plot(time, sratemax_ns, 'dodgerblue', linewidth=2, linestyle='--')
# ax[2].plot(time, pratemax_ss, 'crimson', linewidth=2)
# ax[2].plot(time, sratemax_ss, 'crimson', linewidth=2, linestyle='--')
ax[2].set_xlim([0,32400])
ax[2].set_ylim([0,0.1]) #[20,50] for 5 km
ax[2].set_xlabel('Time (s)', fontsize=14)
ax[2].set_ylabel("Rain rate (kg m$^{-2}$ s$^{-1}$)", fontsize=14)
ax[2].tick_params(axis='both', labelsize=12)
ax[2].set_title(f"Max. surface rain/hail rate", fontsize=16)
ax[2].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[2].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[2].xaxis.set_major_locator(MultipleLocator(3600))
ax[2].xaxis.set_minor_locator(MultipleLocator(900))
ax[2].yaxis.set_major_locator(MultipleLocator(0.02))
ax[2].yaxis.set_minor_locator(MultipleLocator(0.01))

ax3 = ax[2].twinx()
l8,= ax3.plot(time[:-2], movmean(sratemax_fs,5)[:-2], 'k', linewidth=1.5, linestyle='--')
ax3.plot(time[:-2], movmean(sratemax_ns,5)[:-2], 'dodgerblue', linewidth=1.5, linestyle='--')
ax3.plot(time[:-2], movmean(sratemax_ss,5)[:-2], 'crimson', linewidth=1.5, linestyle='--')
ax3.set_ylim([0,0.001])
ax3.set_ylabel("Hail rate (kg m$^{-2}$ s$^{-1}$)", fontsize=14)
ax3.yaxis.set_major_locator(MultipleLocator(0.0002))
ax3.yaxis.set_minor_locator(MultipleLocator(0.0001))
ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax3.tick_params(axis='y', labelsize=12)
# ax[2].legend(handles=[l7,l8,l9], labels=['FREESLIP','NOSLIP','SEMISLIP'],
#              loc='upper left', fontsize=14)
ax[2].legend(handles=[l7,l8], labels=['Rain','Hail'],
             loc='upper left', fontsize=14)

if figsave:
    plt.savefig(fp2+'figs/coldpool_all_timeseries.png', dpi=300)







#%% Time block statistical significance plot


times_hourly = np.linspace(0.5, 8.5, 9)
s = 60

fig,ax = plt.subplots(1, 1, figsize=(12,10.5), layout='constrained')

# wmax1000
ax.plot(np.linspace(0,9,10), 60*np.ones(shape=(10,)), 'red',
        np.linspace(0,9,10), 59*np.ones(shape=(10,)), 'gold',
        np.linspace(0,9,10), 58*np.ones(shape=(10,)), 'mediumblue',
        np.linspace(0,9,10), 57*np.ones(shape=(10,)), 'k', linewidth=1.5)
ax.scatter(times_hourly, 60*np.ma.masked_array(isf['wmax1000'], isf['wmax1000']==0), s=s, marker='o', c='red', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 59*np.ma.masked_array(inf['wmax1000'], inf['wmax1000']==0), s=s, marker='o', c='gold', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 58*np.ma.masked_array(isn['wmax1000'], isn['wmax1000']==0), s=s, marker='o', c='mediumblue', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 57*np.ma.masked_array(ian['wmax1000'], ian['wmax1000']==0), s=s, marker='o', c='k', edgecolors='k', linewidths=1)
# wmax2500
ax.plot(np.linspace(0,9,10), 54*np.ones(shape=(10,)), 'red',
        np.linspace(0,9,10), 53*np.ones(shape=(10,)), 'gold',
        np.linspace(0,9,10), 52*np.ones(shape=(10,)), 'mediumblue',
        np.linspace(0,9,10), 51*np.ones(shape=(10,)), 'k', linewidth=1.5)
ax.scatter(times_hourly, 54*np.ma.masked_array(isf['wmax2500'], isf['wmax2500']==0), s=s, marker='o', c='red', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 53*np.ma.masked_array(inf['wmax2500'], inf['wmax2500']==0), s=s, marker='o', c='gold', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 52*np.ma.masked_array(isn['wmax2500'], isn['wmax2500']==0), s=s, marker='o', c='mediumblue', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 51*np.ma.masked_array(ian['wmax2500'], ian['wmax2500']==0), s=s, marker='o', c='k', edgecolors='k', linewidths=1)
# wmax5000
ax.plot(np.linspace(0,9,10), 48*np.ones(shape=(10,)), 'red',
        np.linspace(0,9,10), 47*np.ones(shape=(10,)), 'gold',
        np.linspace(0,9,10), 46*np.ones(shape=(10,)), 'mediumblue',
        np.linspace(0,9,10), 45*np.ones(shape=(10,)), 'k', linewidth=1.5)
ax.scatter(times_hourly, 48*np.ma.masked_array(isf['wmax5000'], isf['wmax5000']==0), s=s, marker='o', c='red', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 47*np.ma.masked_array(inf['wmax5000'], inf['wmax5000']==0), s=s, marker='o', c='gold', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 46*np.ma.masked_array(isn['wmax5000'], isn['wmax5000']==0), s=s, marker='o', c='mediumblue', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 45*np.ma.masked_array(ian['wmax5000'], ian['wmax5000']==0), s=s, marker='o', c='k', edgecolors='k', linewidths=1)
# vortsfc
ax.plot(np.linspace(0,9,10), 42*np.ones(shape=(10,)), 'red',
        np.linspace(0,9,10), 41*np.ones(shape=(10,)), 'gold',
        np.linspace(0,9,10), 40*np.ones(shape=(10,)), 'mediumblue',
        np.linspace(0,9,10), 39*np.ones(shape=(10,)), 'k', linewidth=1.5)
ax.scatter(times_hourly, 42*np.ma.masked_array(isf['vortsfc'], isf['vortsfc']==0), s=s, marker='o', c='red', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 41*np.ma.masked_array(inf['vortsfc'], inf['vortsfc']==0), s=s, marker='o', c='gold', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 40*np.ma.masked_array(isn['vortsfc'], isn['vortsfc']==0), s=s, marker='o', c='mediumblue', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 39*np.ma.masked_array(ian['vortsfc'], ian['vortsfc']==0), s=s, marker='o', c='k', edgecolors='k', linewidths=1)
# vort1km
ax.plot(np.linspace(0,9,10), 36*np.ones(shape=(10,)), 'red',
        np.linspace(0,9,10), 35*np.ones(shape=(10,)), 'gold',
        np.linspace(0,9,10), 34*np.ones(shape=(10,)), 'mediumblue',
        np.linspace(0,9,10), 33*np.ones(shape=(10,)), 'k', linewidth=1.5)
ax.scatter(times_hourly, 36*np.ma.masked_array(isf['vort1km'], isf['vort1km']==0), s=s, marker='o', c='red', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 35*np.ma.masked_array(inf['vort1km'], inf['vort1km']==0), s=s, marker='o', c='gold', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 34*np.ma.masked_array(isn['vort1km'], isn['vort1km']==0), s=s, marker='o', c='mediumblue', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 33*np.ma.masked_array(ian['vort1km'], ian['vort1km']==0), s=s, marker='o', c='k', edgecolors='k', linewidths=1)
# vort3km
ax.plot(np.linspace(0,9,10), 30*np.ones(shape=(10,)), 'red',
        np.linspace(0,9,10), 29*np.ones(shape=(10,)), 'gold',
        np.linspace(0,9,10), 28*np.ones(shape=(10,)), 'mediumblue',
        np.linspace(0,9,10), 27*np.ones(shape=(10,)), 'k', linewidth=1.5)
ax.scatter(times_hourly, 30*np.ma.masked_array(isf['vort3km'], isf['vort3km']==0), s=s, marker='o', c='red', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 29*np.ma.masked_array(inf['vort3km'], inf['vort3km']==0), s=s, marker='o', c='gold', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 28*np.ma.masked_array(isn['vort3km'], isn['vort3km']==0), s=s, marker='o', c='mediumblue', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 27*np.ma.masked_array(ian['vort3km'], ian['vort3km']==0), s=s, marker='o', c='k', edgecolors='k', linewidths=1)
# swspmax
ax.plot(np.linspace(0,9,10), 24*np.ones(shape=(10,)), 'red',
        np.linspace(0,9,10), 23*np.ones(shape=(10,)), 'gold',
        np.linspace(0,9,10), 22*np.ones(shape=(10,)), 'mediumblue',
        np.linspace(0,9,10), 21*np.ones(shape=(10,)), 'k', linewidth=1.5)
ax.scatter(times_hourly, 24*np.ma.masked_array(isf['swspmax'], isf['swspmax']==0), s=s, marker='o', c='red', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 23*np.ma.masked_array(inf['swspmax'], inf['swspmax']==0), s=s, marker='o', c='gold', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 22*np.ma.masked_array(isn['swspmax'], isn['swspmax']==0), s=s, marker='o', c='mediumblue', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 21*np.ma.masked_array(ian['swspmax'], ian['swspmax']==0), s=s, marker='o', c='k', edgecolors='k', linewidths=1)
# sthpmin
ax.plot(np.linspace(0,9,10), 18*np.ones(shape=(10,)), 'red',
        np.linspace(0,9,10), 17*np.ones(shape=(10,)), 'gold',
        np.linspace(0,9,10), 16*np.ones(shape=(10,)), 'mediumblue',
        np.linspace(0,9,10), 15*np.ones(shape=(10,)), 'k', linewidth=1.5)
ax.scatter(times_hourly, 18*np.ma.masked_array(isf['sthpmin'], isf['sthpmin']==0), s=s, marker='o', c='red', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 17*np.ma.masked_array(inf['sthpmin'], inf['sthpmin']==0), s=s, marker='o', c='gold', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 16*np.ma.masked_array(isn['sthpmin'], isn['sthpmin']==0), s=s, marker='o', c='mediumblue', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 15*np.ma.masked_array(ian['sthpmin'], ian['sthpmin']==0), s=s, marker='o', c='k', edgecolors='k', linewidths=1)
# pratemax
ax.plot(np.linspace(0,9,10), 12*np.ones(shape=(10,)), 'red',
        np.linspace(0,9,10), 11*np.ones(shape=(10,)), 'gold',
        np.linspace(0,9,10), 10*np.ones(shape=(10,)), 'mediumblue',
        np.linspace(0,9,10), 9*np.ones(shape=(10,)), 'k', linewidth=1.5)
ax.scatter(times_hourly, 12*np.ma.masked_array(isf['pratemax'], isf['pratemax']==0), s=s, marker='o', c='red', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 11*np.ma.masked_array(inf['pratemax'], inf['pratemax']==0), s=s, marker='o', c='gold', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 10*np.ma.masked_array(isn['pratemax'], isn['pratemax']==0), s=s, marker='o', c='mediumblue', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 9*np.ma.masked_array(ian['pratemax'], ian['pratemax']==0), s=s, marker='o', c='k', edgecolors='k', linewidths=1)
# sratemax
l1,= ax.plot(np.linspace(0,9,10), 6*np.ones(shape=(10,)), 'red', linewidth=1.5)
l2,= ax.plot(np.linspace(0,9,10), 5*np.ones(shape=(10,)), 'gold', linewidth=1.5)
l3,= ax.plot(np.linspace(0,9,10), 4*np.ones(shape=(10,)), 'mediumblue', linewidth=1.5)
l4,= ax.plot(np.linspace(0,9,10), 3*np.ones(shape=(10,)), 'k', linewidth=1.5)
ax.scatter(times_hourly, 6*np.ma.masked_array(isf['sratemax'], isf['sratemax']==0), s=s, marker='o', c='red', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 5*np.ma.masked_array(inf['sratemax'], inf['sratemax']==0), s=s, marker='o', c='gold', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 4*np.ma.masked_array(isn['sratemax'], isn['sratemax']==0), s=s, marker='o', c='mediumblue', edgecolors='k', linewidths=1)
ax.scatter(times_hourly, 3*np.ma.masked_array(ian['sratemax'], ian['sratemax']==0), s=s, marker='o', c='k', edgecolors='k', linewidths=1)
s1 = ax.scatter([-1], [-1], s=50, marker='o', c='w', edgecolors='k', linewidths=1)

ax.grid(visible=True, which='minor', axis='x', color='darkgray', linestyle='-')
ax.xaxis.set_minor_locator(MultipleLocator(1))
xlab = ['0-1 h', '1-2 h', '2-3 h', '3-4 h', '4-5 h', '5-6 h', '6-7 h', '7-8 h', '8-9 h']
ylab = ["Hail\n rate", "Rain\n rate", "10-m \u03B8'", "10-m WS", "3-km \u03B6", "1-km \u03B6", "10-m \u03B6", "5-km w", "2.5-km w", "1-km w"]
ax.set_xlim([0,9])
ax.set_ylim([0,68])
ax.set_xticks(ticks=times_hourly, labels=xlab, fontsize=14)
ax.set_yticks(ticks=np.linspace(4.5, 58.5, 10), labels=ylab, fontsize=14)
ax.set_xlabel('Time (h)', fontsize=16)
ax.legend(handles=[l1,l2,l3,l4,s1], labels=['SS-FS','NS-FS','SS-NS','ANOVA','p<0.05'], ncol=2,
             loc='upper right', fontsize=13)
ax.text(0.1, 65.5, "Hourly statistical significance", fontsize=25, fontweight='bold')



figsave = False

if figsave:
    plt.savefig(fp2+"figs/stat_sig.png", dpi=300)






#%% Swaths?

# Freeslip
fp1 = 'C:/Users/mschne28/Documents/cm1out/cwe/freeslip_wk_250m/'
ds = nc.Dataset(fp1+'cm1out_000037.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
sws_fs = ds.variables['sws'][:].data[0,:,:] #max sfc wind
svs_fs = ds.variables['svs'][:].data[0,:,:] #max sfc vort
sus_fs = ds.variables['sus'][:].data[0,:,:] #max 5km updraft
shs_fs = ds.variables['shs'][:].data[0,:,:] #max integrated UH
ds.close()


# Semislip
fp2 = 'C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_250m/'
ds = nc.Dataset(fp2+'cm1out_000037.nc')
sws_ss = ds.variables['sws'][:].data[0,:,:]
svs_ss = ds.variables['svs'][:].data[0,:,:]
sus_ss = ds.variables['sus'][:].data[0,:,:]
shs_ss = ds.variables['shs'][:].data[0,:,:]
ds.close()

# Noslip
fp3 = 'C:/Users/mschne28/Documents/cm1out/cwe/noslip_wk_250m/'
ds = nc.Dataset(fp3+'cm1out_000037.nc')
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


#%% Assemble translated swaths



fp1 = 'C:/Users/mschne28/Documents/cm1out/cwe/freeslip_wk_250m/'
fp2 = 'C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_250m/'
fp3 = 'C:/Users/mschne28/Documents/cm1out/cwe/noslip_wk_250m/'

fn = np.linspace(5,37,9)


ds = nc.Dataset(fp1+"cm1out_000001.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
umove = ds.variables['umove'][:].data
vmove = ds.variables['vmove'][:].data
ds.close()

# umove = 20
# vmove = 2


x_added = umove*3600/1000
nx_added = np.round(x_added/0.25)
nxt_added = nx_added*(len(fn)-1)
xt_added = nxt_added*0.25
xn = np.arange(xh[0], xh[-1]+xt_added+0.25, 0.25)


y_added = vmove*3600/1000 # y distance added per hour
ny_added = np.round(y_added/0.25) # number of y grid points added per hour
nyt_added = ny_added*(len(fn)-1) # total number y grid points added over 8 h
yt_added = nyt_added*0.25  # total y distance added over 8 h (km)
yn = np.arange(yh[0], yh[-1]+yt_added+0.25, 0.25)


sws_fs = np.zeros(shape=(len(yn), len(xn)), dtype=float); svs_fs = np.zeros(shape=(len(yn), len(xn)), dtype=float)
sus_fs = np.zeros(shape=(len(yn), len(xn)), dtype=float); shs_fs = np.zeros(shape=(len(yn), len(xn)), dtype=float)
sws_ss = np.zeros(shape=(len(yn), len(xn)), dtype=float); svs_ss = np.zeros(shape=(len(yn), len(xn)), dtype=float)
sus_ss = np.zeros(shape=(len(yn), len(xn)), dtype=float); shs_ss = np.zeros(shape=(len(yn), len(xn)), dtype=float)
sws_ns = np.zeros(shape=(len(yn), len(xn)), dtype=float); svs_ns = np.zeros(shape=(len(yn), len(xn)), dtype=float)
sus_ns = np.zeros(shape=(len(yn), len(xn)), dtype=float); shs_ns = np.zeros(shape=(len(yn), len(xn)), dtype=float)





for f in fn:
    n = (f-fn[0])/(fn[1]-fn[0])
    
    ix = slice(int(nx_added*n), int(nx_added*n + len(xh)))
    iy = slice(int(ny_added*n), int(ny_added*n + len(yh)))
    
    # Freeslip
    ds = nc.Dataset(fp1+f"cm1out_{f:06.0f}.nc")
    sws2 = ds.variables['sws2'][:].data[0,:,:] #max sfc wind
    svs2 = ds.variables['svs2'][:].data[0,:,:] #max sfc vort
    sus2 = ds.variables['sus2'][:].data[0,:,:] #max 5km updraft
    shs2 = ds.variables['shs2'][:].data[0,:,:] #max integrated UH
    
    sws_fs[iy,ix] = np.maximum(sws_fs[iy,ix], sws2)
    svs_fs[iy,ix] = np.maximum(svs_fs[iy,ix], svs2)
    sus_fs[iy,ix] = np.maximum(sus_fs[iy,ix], sus2)
    shs_fs[iy,ix] = np.maximum(shs_fs[iy,ix], shs2)
    ds.close()
    
    
    # Semislip
    ds = nc.Dataset(fp2+f"cm1out_{f:06.0f}.nc")
    sws2 = ds.variables['sws2'][:].data[0,:,:]
    svs2 = ds.variables['svs2'][:].data[0,:,:]
    sus2 = ds.variables['sus2'][:].data[0,:,:]
    shs2 = ds.variables['shs2'][:].data[0,:,:]
    
    sws_ss[iy,ix] = np.maximum(sws_ss[iy,ix], sws2)
    svs_ss[iy,ix] = np.maximum(svs_ss[iy,ix], svs2)
    sus_ss[iy,ix] = np.maximum(sus_ss[iy,ix], sus2)
    shs_ss[iy,ix] = np.maximum(shs_ss[iy,ix], shs2)
    ds.close()
    
    
    # Noslip
    ds = nc.Dataset(fp3+f"cm1out_{f:06.0f}.nc")
    sws2 = ds.variables['sws2'][:].data[0,:,:]
    svs2 = ds.variables['svs2'][:].data[0,:,:]
    sus2 = ds.variables['sus2'][:].data[0,:,:]
    shs2 = ds.variables['shs2'][:].data[0,:,:]
    
    sws_ns[iy,ix] = np.maximum(sws_ns[iy,ix], sws2)
    svs_ns[iy,ix] = np.maximum(svs_ns[iy,ix], svs2)
    sus_ns[iy,ix] = np.maximum(sus_ns[iy,ix], sus2)
    shs_ns[iy,ix] = np.maximum(shs_ns[iy,ix], shs2)
    ds.close()
    

# if False:
#     dbfile = open(fp1+"composite_swaths_30min.pkl", 'wb')
#     data1 = {'xn':xn, 'yn':yn, 'sws':sws_fs, 'svs':svs_fs, 'sus':sus_fs, 'shs':shs_fs}
#     pickle.dump(data1, dbfile)
#     dbfile.close()
        
#     dbfile = open(fp2+"composite_swaths_30min.pkl", 'wb')
#     data2 = {'xn':xn, 'yn':yn, 'sws':sws_ss, 'svs':svs_ss, 'sus':sus_ss, 'shs':shs_ss}
#     pickle.dump(data2, dbfile)
#     dbfile.close()
    
#     dbfile = open(fp3+"composite_swaths_30min.pkl", 'wb')
#     data3 = {'xn':xn, 'yn':yn, 'sws':sws_ns, 'svs':svs_ns, 'sus':sus_ns, 'shs':shs_ns}
#     pickle.dump(data3, dbfile)
#     dbfile.close()

    
#%% Plot translated swaths


dbfile = open(fp1+"composite_swaths_30min.pkl", 'rb')
fs = pickle.load(dbfile)
xn = fs['xn']
yn = fs['yn']
sws_fs = fs['sws']
svs_fs = fs['svs']
sus_fs = fs['sus']
shs_fs = fs['shs']
dbfile.close()

dbfile = open(fp2+"composite_swaths_30min.pkl", 'rb')
ss = pickle.load(dbfile)
sws_ss = ss['sws']
svs_ss = ss['svs']
sus_ss = ss['sus']
shs_ss = ss['shs']
dbfile.close()

dbfile = open(fp3+"composite_swaths_30min.pkl", 'rb')
ns = pickle.load(dbfile)
sws_ns = ns['sws']
svs_ns = ns['svs']
sus_ns = ns['sus']
shs_ns = ns['shs']
dbfile.close()

xn2 = xn + 150
yn2 = yn + 150


xl = [0,900]
yl = [0,360]


figsave = False


# 5-km updraft swaths
fig,ax = plt.subplots(3, 1, figsize=(10,11), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(aspect='equal'))

plot_contourf(xn2, yn2, sus_fs, 'w', ax[0],
              levels=np.linspace(0,40,21), datalims=[0,40],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max')
# ax[0].contour(xn2, yn2, svs_fs, levels=[0.04], colors='k', linewidths=[0.75])
# l1, = ax[0].plot([-190,-200], [-190,-200], '-k', linewidth=0.75)
# ax[0].legend(handles=[l1], labels=["\u03B6=0.04 s$^{-1}$"], loc='lower left', fontsize=12)
ax[0].set_title("5-km updraft", fontsize=16)
ax[0].set_ylabel('y (km)', fontsize=14)


plot_contourf(xn2, yn2, sus_ss, 'w', ax[1],
              levels=np.linspace(0,40,21), datalims=[0,40],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max')
# ax[1].contour(xn2, yn2, svs_ss, levels=[0.04], colors='k', linewidths=[0.75])
ax[1].set_ylabel('y (km)', fontsize=14)


plot_contourf(xn2, yn2, sus_ns, 'w', ax[2],
              levels=np.linspace(0,40,21), datalims=[0,40],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max')
# ax[2].contour(xn2, yn2, svs_ns, levels=[0.04], colors='k', linewidths=[0.75])
ax[2].set_xlabel('x (km)', fontsize=14)
ax[2].set_ylabel('y (km)', fontsize=14)

if figsave:
    plt.savefig(fp2+'figs/w5km_swaths.png', dpi=300)



# Sfc wind speed swaths
fig,ax = plt.subplots(3, 1, figsize=(10,11), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(aspect='equal'))

plot_contourf(xn2, yn2, sws_fs, 'wspd', ax[0],
              levels=np.linspace(0,40,21), datalims=[0,40],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max')
ax[0].contour(xn2, yn2, sws_fs, levels=[26,33], colors=['navy','k'], linewidths=[0.5,1])
l1, = ax[0].plot([-190,-200], [-190,-200], color='navy', linewidth=0.5)
l2, = ax[0].plot([-190,-200], [-190,-200], '-k', linewidth=1)
ax[0].legend(handles=[l1,l2], labels=['Severe (26 m/s)','Sig. Severe (33 m/s)'], loc='lower left', fontsize=12)
ax[0].set_title("Surface wind speed", fontsize=16)
ax[0].set_ylabel('y (km)', fontsize=14)


plot_contourf(xn2, yn2, sws_ss, 'wspd', ax[1],
              levels=np.linspace(0,40,21), datalims=[0,40],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max')
ax[1].contour(xn2, yn2, sws_ss, levels=[26,33], colors=['navy','k'], linewidths=[0.5,1])
ax[1].set_ylabel('y (km)', fontsize=14)


plot_contourf(xn2, yn2, sws_ns, 'wspd', ax[2],
              levels=np.linspace(0,40,21), datalims=[0,40],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max')
ax[2].contour(xn2, yn2, sws_ns, levels=[26,33], colors=['navy','k'], linewidths=[0.5,1])
ax[2].set_xlabel('x (km)', fontsize=14)
ax[2].set_ylabel('y (km)', fontsize=14)

if figsave:
    plt.savefig(fp2+'figs/wspd_swaths.png', dpi=300)



# Sfc vorticity swaths
fig,ax = plt.subplots(3, 1, figsize=(10,11), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(aspect='equal'))

plot_contourf(xn2, yn2, svs_fs, 'zvort', ax[0],
              levels=np.linspace(0,0.05,21), datalims=[0,0.05],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,0.05,11), extend='max')
# ax[0].contour(xn2, yn2, sus_fs, levels=[20,35], colors=['dimgray','k'], linewidths=[0.75,1.5])
# l2, = ax[0].plot([-190,-200], [-190,-200], color='dimgray', linewidth=0.75)
# l3, = ax[0].plot([-190,-200], [-190,-200], '-k', linewidth=1.5)
# ax[0].legend(handles=[l2,l3], labels=['w=20 m/s','w=35 m/s'], loc='upper right', fontsize=12)
ax[0].contour(xn2, yn2, sws_fs, levels=[26,33], colors=['navy','k'], linewidths=[0.5,1])
l3, = ax[0].plot([-190,-200], [-190,-200], color='navy', linewidth=0.5)
l4, = ax[0].plot([-190,-200], [-190,-200], '-k', linewidth=1)
ax[0].legend(handles=[l3,l4], labels=['Severe (26 m/s)','Sig. Severe (33 m/s)'], loc='lower left', fontsize=12)
ax[0].set_title("Surface vorticity", fontsize=16)
ax[0].set_ylabel('y (km)', fontsize=14)

plot_contourf(xn2, yn2, svs_ss, 'zvort', ax[1],
              levels=np.linspace(0,0.05,21), datalims=[0,0.05],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,0.05,11), extend='max')
# ax[1].contour(xn2, yn2, sus_ss, levels=[20,35], colors=['dimgray','k'], linewidths=[0.75,1.5])
ax[1].contour(xn2, yn2, sws_ss, levels=[26,33], colors=['navy','k'], linewidths=[0.5,1])
ax[1].set_ylabel('y (km)', fontsize=14)

plot_contourf(xn2, yn2, svs_ns, 'zvort', ax[2],
              levels=np.linspace(0,0.05,21), datalims=[0,0.05],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,0.05,11), extend='max')
# ax[2].contour(xn2, yn2, sus_ns, levels=[20,35], colors=['dimgray','k'], linewidths=[0.75,1.5])
ax[2].contour(xn2, yn2, sws_ns, levels=[26,33], colors=['navy','k'], linewidths=[0.5,1])
ax[2].set_xlabel('x (km)', fontsize=14)
ax[2].set_ylabel('y (km)', fontsize=14)

if figsave:
    plt.savefig(fp2+'figs/vort_swaths.png', dpi=300)



# Max updraft helicity swaths
fig,ax = plt.subplots(3, 1, figsize=(10,11), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(aspect='equal'))

plot_contourf(xn2, yn2, shs_fs, 'uh', ax[0],
              levels=np.linspace(0,1600,17), datalims=[0,1600],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,1600,9), extend='max')
# ax[0].contour(xn2, yn2, sus_fs, levels=[25,35], colors=['dimgray','k'], linewidths=[0.75,1.5])
# l5, = ax[0].plot([-190,-200], [-190,-200], color='dimgray', linewidth=0.75)
# l6, = ax[0].plot([-190,-200], [-190,-200], '-k', linewidth=1.5)
# ax[0].legend(handles=[l5,l6], labels=['w=25 m/s','w=35 m/s'], loc='lower left', fontsize=12)
ax[0].set_title("Updraft helicity", fontsize=16)
ax[0].set_ylabel('y (km)', fontsize=14)

plot_contourf(xn2, yn2, shs_ss, 'uh', ax[1],
              levels=np.linspace(0,1600,17), datalims=[0,1600],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,1600,9), extend='max')
# ax[1].contour(xn2, yn2, sus_ss, levels=[25,35], colors=['dimgray','k'], linewidths=[0.75,1.5])
ax[1].set_ylabel('y (km)', fontsize=14)

plot_contourf(xn2, yn2, shs_ns, 'uh', ax[2],
              levels=np.linspace(0,1600,17), datalims=[0,1600],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,1600,9), extend='max')
# ax[2].contour(xn2, yn2, sus_ns, levels=[25,35], colors=['dimgray','k'], linewidths=[0.75,1.5])
ax[2].set_xlabel('x (km)', fontsize=14)
ax[2].set_ylabel('y (km)', fontsize=14)

if figsave:
    plt.savefig(fp2+'figs/UH_swaths.png', dpi=300)


plt.show()



#%% Comparison between P3 and NSSL (and Morrison if it ever runs)

fp1 = 'C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_500m/'
fp2 = 'C:/Users/mschne28/Documents/cm1out/cwe/semislip_nssl_500m/'
# fp3 = 'C:/Users/mschne28/Documents/cm1out/cwe/semislip_morr_500m/'

fn = np.linspace(5,37,5)


figsave = False

# fig1,ax1 = plt.subplots(2, 5, figsize=(12.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
# fig2,ax2 = plt.subplots(2, 5, figsize=(12.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
fig1,ax1 = plt.subplots(2, 5, figsize=(12.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
fig2,ax2 = plt.subplots(2, 5, figsize=(12.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')


for f in fn:
    ### P3 scheme
    ds = nc.Dataset(fp1+f"cm1out_{f:06.0f}.nc")
    time = ds.variables['time'][:].data[0]
    xh = ds.variables['xh'][:].data
    yh = ds.variables['yh'][:].data
    zh = ds.variables['zh'][:].data
    iz1 = np.where(zh>1)[0][0]
    iz2 = np.where(zh>2)[0][0]
    iz3 = np.where(zh>3)[0][0]
    
    
    dbz = ds.variables['dbz'][:].data[0,0,:,:]
    winterp = ds.variables['winterp'][:].data[0,0:iz2,:,:]
    zvort = ds.variables['zvort'][:].data[0,iz1:iz2,:,:]
    thrpert = ds.variables['th'][:].data[0,0,:,:] - ds.variables['th0'][:].data[0,0,:,:]
    uinterp = ds.variables['uinterp'][:].data[0,0,:,:]
    vinterp = ds.variables['vinterp'][:].data[0,0,:,:]
    u_gr = uinterp + ds.variables['umove'][:].data[0]
    v_gr = vinterp + ds.variables['vmove'][:].data[0]
    # ### P3 3-moment scheme
    # if 'qi1' in list(ds.variables.keys()):
    #     thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] - 
    #                 (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] + 
    #                  ds.variables['qi1'][:].data[0,0,:,:] +
    #                  ds.variables['qi2'][:].data[0,0,:,:] + 
    #                  ds.variables['qi3'][:].data[0,0,:,:]))
    #                  # ds.variables['qi4'][:].data[0,0,:,:]))
    # ### NSSL 3-moment scheme
    # elif 'qg' in list(ds.variables.keys()):
    #     thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] -
    #                 (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] +
    #                  ds.variables['qi'][:].data[0,0,:,:] + ds.variables['qs'][:].data[0,0,:,:] +
    #                  ds.variables['qg'][:].data[0,0,:,:] + ds.variables['qhl'][:].data[0,0,:,:]))
    
    # thr0 = ds.variables['th0'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,:,:])
    # thpert = thr - thr0
    # del thr,thr0
    ds.close()
    
    
    xl = [-150,150]
    yl = [-150,150]
    
    # xl = [-100,100]
    # yl = [-100,100]
    
    
    n = int((f-fn[0])/(fn[1]-fn[0]))
    
    if f == fn[-1]:
        cb_flag = True
    else:
        cb_flag = False
    
    
    qix = 30
    
    
    plot_contourf(xh, yh, np.ma.masked_array(dbz, dbz<0.1), 'dbz', ax1[0,n], levels=np.linspace(0,70,15),
                  datalims=[0,70], xlims=xl, ylims=yl, cmap='HomeyerRainbow', cbar=cb_flag, cbfs=10)
    ax1[0,n].contour(xh, yh, np.max(winterp, axis=0), levels=[5,10], colors=['dimgray','k'], linestyles='-', linewidths=[0.75,0.75])
    if n == 0:
        l1, = ax1[0,0].plot([190,200], [190,200], color='dimgray', linewidth=0.75)
        l2, = ax1[0,0].plot([190,200], [190,200], '-k', linewidth=0.75)
        ax1[0,0].legend(handles=[l1,l2], labels=['w=5 m/s','w=10 m/s'], loc='upper right', fontsize=10)
    ax1[0,n].set_title(f"t = {time:.0f} s")
    # fig1.suptitle(f"Sfc dbz + max 0-2 km w ({titlestr})")
    
    
    
    plot_contourf(xh, yh, thrpert, 'thpert', ax2[0,n], levels=np.linspace(-12,12,25),
                  datalims=[-12,12], xlims=xl, ylims=yl, cmap='balance', cbar=cb_flag, cbfs=10)
    ax2[0,n].contour(xh, yh, np.max(zvort, axis=0), levels=[0.015], colors='r', linestyles='-', linewidths=1)
    ax2[0,n].quiver(xh[::qix], yh[::qix], u_gr[::qix,::qix], v_gr[::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
    if n == 0:
        l3, = ax2[0,0].plot([190,200], [190,200], '-r', linewidth=1)
        ax2[0,0].legend(handles=[l3], labels=["\u03B6=0.015 s$^{-1}$"], loc='upper right', fontsize=10)
    ax2[0,n].set_title(f"t = {time:.0f} s")
    # fig2.suptitle(f"Sfc thrpert + sfc wind + max 0-1 km zeta=0.025 s$^{{-1}}$ ({titlestr})")
    
    
    
    ### NSSL scheme
    ds = nc.Dataset(fp2+f"cm1out_{f:06.0f}.nc")
    dbz = ds.variables['dbz'][:].data[0,0,:,:]
    winterp = ds.variables['winterp'][:].data[0,0:iz2,:,:]
    zvort = ds.variables['zvort'][:].data[0,iz1:iz2,:,:]
    thrpert = ds.variables['th'][:].data[0,0,:,:] - ds.variables['th0'][:].data[0,0,:,:]
    uinterp = ds.variables['uinterp'][:].data[0,0,:,:]
    vinterp = ds.variables['vinterp'][:].data[0,0,:,:]
    u_gr = uinterp + ds.variables['umove'][:].data[0]
    v_gr = vinterp + ds.variables['vmove'][:].data[0]
    # ### P3 3-moment scheme
    # if 'qi1' in list(ds.variables.keys()):
    #     thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] - 
    #                 (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] + 
    #                  ds.variables['qi1'][:].data[0,0,:,:] +
    #                  ds.variables['qi2'][:].data[0,0,:,:] + 
    #                  ds.variables['qi3'][:].data[0,0,:,:]))
    #                  # ds.variables['qi4'][:].data[0,0,:,:]))
    # ### NSSL 3-moment scheme
    # if 'qg' in list(ds.variables.keys()):
    #     thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] -
    #                 (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] +
    #                  ds.variables['qi'][:].data[0,0,:,:] + ds.variables['qs'][:].data[0,0,:,:] +
    #                  ds.variables['qg'][:].data[0,0,:,:] + ds.variables['qhl'][:].data[0,0,:,:]))
    
    # thr0 = ds.variables['th0'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,:,:])
    # thpert = thr - thr0
    # del thr,thr0
    ds.close()
    
    
    plot_contourf(xh, yh, np.ma.masked_array(dbz, dbz<0.1), 'dbz', ax1[1,n], levels=np.linspace(0,70,15),
                  datalims=[0,70], xlims=xl, ylims=yl, cmap='HomeyerRainbow', cbar=cb_flag, cbfs=10)
    ax1[1,n].contour(xh, yh, np.max(winterp, axis=0), levels=[5,10], colors=['dimgray','k'], linestyles='-', linewidths=[0.75,0.75])
    # ax1[1,n].set_title(f"t = {time:.0f} s")
    # fig1.suptitle(f"Sfc dbz + max 0-2 km w ({titlestr})")
    # if (n==len(fn)-1) & (figsave):
    #     fig1.savefig(fp1+f"figs/dbz_compare.png", dpi=300)
    
    
    
    plot_contourf(xh, yh, thrpert, 'thpert', ax2[1,n], levels=np.linspace(-12,12,25),
                  datalims=[-12,12], xlims=xl, ylims=yl, cmap='balance', cbar=cb_flag, cbfs=10)
    ax2[1,n].contour(xh, yh, np.max(zvort, axis=0), levels=[0.015], colors='r', linestyles='-', linewidths=1)
    ax2[1,n].quiver(xh[::qix], yh[::qix], u_gr[::qix,::qix], v_gr[::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
    # ax2[1,n].set_title(f"t = {time:.0f} s")
    # fig2.suptitle(f"Sfc thrpert + sfc wind + max 0-1 km zeta=0.025 s$^{{-1}}$ ({titlestr})")
    # if (n==len(fn)-1) & (figsave):
    #     fig2.savefig(fp1+f"figs/thrpert_compare.png", dpi=300)
    
    
    
    
    ### Morrison scheme
    ds = nc.Dataset(fp3+f"cm1out_{f:06.0f}.nc")
    dbz = ds.variables['dbz'][:].data[0,0,:,:]
    winterp = ds.variables['winterp'][:].data[0,0:iz2,:,:]
    zvort = ds.variables['zvort'][:].data[0,iz1:iz2,:,:]
    thrpert = ds.variables['th'][:].data[0,0,:,:] - ds.variables['th0'][:].data[0,0,:,:]
    uinterp = ds.variables['uinterp'][:].data[0,0,:,:]
    vinterp = ds.variables['vinterp'][:].data[0,0,:,:]
    u_gr = uinterp + ds.variables['umove'][:].data[0]
    v_gr = vinterp + ds.variables['vmove'][:].data[0]
    # ### P3 3-moment scheme
    # if 'qi1' in list(ds.variables.keys()):
    #     thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] - 
    #                 (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] + 
    #                  ds.variables['qi1'][:].data[0,0,:,:] +
    #                  ds.variables['qi2'][:].data[0,0,:,:] + 
    #                  ds.variables['qi3'][:].data[0,0,:,:]))
    #                  # ds.variables['qi4'][:].data[0,0,:,:]))
    # ### NSSL 3-moment scheme
    # if 'qg' in list(ds.variables.keys()):
    #     thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] -
    #                 (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] +
    #                  ds.variables['qi'][:].data[0,0,:,:] + ds.variables['qs'][:].data[0,0,:,:] +
    #                  ds.variables['qg'][:].data[0,0,:,:] + ds.variables['qhl'][:].data[0,0,:,:]))
    
    # thr0 = ds.variables['th0'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,:,:])
    # thpert = thr - thr0
    # del thr,thr0
    ds.close()
    
    
    plot_contourf(xh, yh, np.ma.masked_array(dbz, dbz<0.1), 'dbz', ax1[2,n], levels=np.linspace(0,70,15),
                  datalims=[0,70], xlims=xl, ylims=yl, cmap='HomeyerRainbow', cbar=cb_flag, cbfs=10)
    ax1[2,n].contour(xh, yh, np.max(winterp, axis=0), levels=[5,10], colors=['dimgray','k'], linestyles='-', linewidths=[0.75,0.75])
    # ax1[2,n].set_title(f"t = {time:.0f} s")
    # fig1.suptitle(f"Sfc dbz + max 0-2 km w ({titlestr})")
    if (n==len(fn)-1) & (figsave):
        fig1.savefig(fp1+f"figs/dbz_compare_v2.png", dpi=300)
    
    
    
    plot_contourf(xh, yh, thrpert, 'thpert', ax2[2,n], levels=np.linspace(-12,12,25),
                  datalims=[-12,12], xlims=xl, ylims=yl, cmap='balance', cbar=cb_flag, cbfs=10)
    ax2[2,n].contour(xh, yh, np.max(zvort, axis=0), levels=[0.015], colors='r', linestyles='-', linewidths=1)
    ax2[2,n].quiver(xh[::qix], yh[::qix], u_gr[::qix,::qix], v_gr[::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
    # ax2[2,n].set_title(f"t = {time:.0f} s")
    # fig2.suptitle(f"Sfc thrpert + sfc wind + max 0-1 km zeta=0.025 s$^{{-1}}$ ({titlestr})")
    if (n==len(fn)-1) & (figsave):
        fig2.savefig(fp1+f"figs/thrpert_compare_v2.png", dpi=300)
    

#%% Get maxima

fp = 'C:/Users/mschne28/Documents/cm1out/cwe/freeslip_wk_250m/'


vortmax_50m = np.zeros(shape=(37,), dtype=float)
vortmax_100m = np.zeros(shape=(37,), dtype=float)
rainmax_sfc = np.zeros(shape=(37,), dtype=float)
hailmax_sfc = np.zeros(shape=(37,), dtype=float)


xy_zv_sfc = np.zeros(shape=(37,2), dtype=float)
xy_zv_50m = np.zeros(shape=(37,2), dtype=float)
xy_zv_100m = np.zeros(shape=(37,2), dtype=float)
xy_zv_1km = np.zeros(shape=(37,2), dtype=float)
xy_zv_3km = np.zeros(shape=(37,2), dtype=float)
xy_w_1km = np.zeros(shape=(37,2), dtype=float)
xy_w_25km = np.zeros(shape=(37,2), dtype=float)
xy_w_5km = np.zeros(shape=(37,2), dtype=float)
xy_wspd_sfc = np.zeros(shape=(37,2), dtype=float)
xy_thp_sfc = np.zeros(shape=(37,2), dtype=float)
xy_rain_sfc = np.zeros(shape=(37,2), dtype=float)
xy_hail_sfc = np.zeros(shape=(37,2), dtype=float)
xy_prate_sfc = np.zeros(shape=(37,2), dtype=float)
xy_srate_sfc = np.zeros(shape=(37,2), dtype=float)



for fn in np.linspace(1,37,37):
    print(f"cm1out_{fn:06.0f}.nc")
    ds = nc.Dataset(fp+f"cm1out_{fn:06.0f}.nc")
    xh = ds.variables['xh'][:].data
    yh = ds.variables['yh'][:].data
    zh = ds.variables['zh'][:].data
    
    
    iz50 = np.argmin(abs(zh-0.05))
    iz90 = np.argmin(abs(zh-0.09))
    iz110 = np.argmin(abs(zh-0.11))
    iz1000 = np.argmin(abs(zh-1.0))
    iz2500 = np.argmin(abs(zh-2.5))
    iz3000 = np.argmin(abs(zh-3.0))
    iz5000 = np.argmin(abs(zh-5.0))
    
    zvort_sfc = ds.variables['zvort'][:].data[0,0,:,:]
    zvort_50m = ds.variables['zvort'][:].data[0,iz50,:,:]
    zvort_100m = np.mean(ds.variables['zvort'][:].data[0,iz90:iz110+1,:,:], axis=0)
    zvort_1000 = ds.variables['zvort'][:].data[0,iz1000,:,:]
    zvort_3000 = ds.variables['zvort'][:].data[0,iz3000,:,:]
    
    w_1000 = ds.variables['winterp'][:].data[0,iz1000,:,:]
    w_2500 = ds.variables['winterp'][:].data[0,iz2500,:,:]
    w_5000 = ds.variables['winterp'][:].data[0,iz5000,:,:]
    
    wspd_sfc = np.sqrt(ds.variables['uinterp'][:].data[0,0,:,:]**2 + ds.variables['vinterp'][:].data[0,0,:,:]**2)
    # wsgr_sfc = np.sqrt((ds.variables['uinterp'][:].data[0,0,:,:]+ds.variables['umove'][:].data)**2 + (ds.variables['vinterp'][:].data[0,0,:,:]+ds.variables['umove'][:].data)**2)
    thp_sfc = ds.variables['th'][:].data[0,0,:,:] - ds.variables['th0'][:].data[0,0,:,:]
    rain_sfc = ds.variables['rain'][:].data[0,:,:]
    hail_sfc = ds.variables['hail'][:].data[0,:,:]
    prate_sfc = ds.variables['prate'][:].data[0,:,:]
    srate_sfc = ds.variables['srate'][:].data[0,:,:]
    
    ds.close()
    
    f = int(fn-1)
    
    vortmax_50m[f] = np.max(zvort_50m)
    vortmax_100m[f] = np.max(zvort_100m)
    rainmax_sfc[f] = np.max(rain_sfc)
    hailmax_sfc[f] = np.max(hail_sfc)
    
    
    i = np.where(zvort_sfc == np.max(zvort_sfc))
    xy_zv_sfc[f,0] = xh[i[1][0]]
    xy_zv_sfc[f,1] = yh[i[0][0]]
    
    i = np.where(zvort_50m == np.max(zvort_50m))
    xy_zv_50m[f,0] = xh[i[1][0]]
    xy_zv_50m[f,1] = yh[i[0][0]]
    
    i = np.where(zvort_100m == np.max(zvort_100m))
    xy_zv_100m[f,0] = xh[i[1][0]]
    xy_zv_100m[f,1] = yh[i[0][0]]
    
    i = np.where(zvort_1000 == np.max(zvort_1000))
    xy_zv_1km[f,0] = xh[i[1][0]]
    xy_zv_1km[f,1] = yh[i[0][0]]
    
    i = np.where(zvort_3000 == np.max(zvort_3000))
    xy_zv_3km[f,0] = xh[i[1][0]]
    xy_zv_3km[f,1] = yh[i[0][0]]
    
    i = np.where(w_1000 == np.max(w_1000))
    xy_w_1km[f,0] = xh[i[1][0]]
    xy_w_1km[f,1] = yh[i[0][0]]
    
    i = np.where(w_2500 == np.max(w_2500))
    xy_w_25km[f,0] = xh[i[1][0]]
    xy_w_25km[f,1] = yh[i[0][0]]
    
    i = np.where(w_5000 == np.max(w_5000))
    xy_w_5km[f,0] = xh[i[1][0]]
    xy_w_5km[f,1] = yh[i[0][0]]
    
    i = np.where(wspd_sfc == np.max(wspd_sfc))
    xy_wspd_sfc[f,0] = xh[i[1][0]]
    xy_wspd_sfc[f,1] = yh[i[0][0]]
    
    i = np.where(thp_sfc == np.min(thp_sfc))
    xy_thp_sfc[f,0] = xh[i[1][0]]
    xy_thp_sfc[f,1] = yh[i[0][0]]
    
    i = np.where(rain_sfc == np.max(rain_sfc))
    xy_rain_sfc[f,0] = xh[i[1][0]]
    xy_rain_sfc[f,1] = yh[i[0][0]]
    
    i = np.where(hail_sfc == np.max(hail_sfc))
    xy_hail_sfc[f,0] = xh[i[1][0]]
    xy_hail_sfc[f,1] = yh[i[0][0]]
    
    i = np.where(prate_sfc == np.max(prate_sfc))
    xy_prate_sfc[f,0] = xh[i[1][0]]
    xy_prate_sfc[f,1] = yh[i[0][0]]
    
    i = np.where(srate_sfc == np.max(srate_sfc))
    xy_srate_sfc[f,0] = xh[i[1][0]]
    xy_srate_sfc[f,1] = yh[i[0][0]]
    
    
    
stat_xy = {'wmax1000':xy_w_1km, 'wmax2500':xy_w_25km, 'wmax5000':xy_w_5km,
           'vortsfc':xy_zv_sfc, 'vort1km':xy_zv_1km, 'vort3km':xy_zv_3km,
           'swpsmax':xy_wspd_sfc, 'sthpmin':xy_thp_sfc, 'pratemax':xy_prate_sfc, 'sratemax':xy_srate_sfc,
           'vort50m':xy_zv_50m, 'vort100m':xy_zv_100m, 'rainmax':xy_rain_sfc, 'hailmax':xy_hail_sfc}


more_stats = {'vort50m':vortmax_50m, 'vort100m':vortmax_100m, 'rainmax':rainmax_sfc, 'hailmax':hailmax_sfc}



dbfile = open(fp+'extra_stats.pkl', 'wb')
pickle.dump(more_stats, dbfile)
dbfile.close()

dbfile = open(fp+'stat_xy.pkl', 'wb')
pickle.dump(stat_xy, dbfile)
dbfile.close()




#%%




















