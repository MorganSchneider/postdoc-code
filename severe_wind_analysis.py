# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 09:57:54 2025

@author: mschne28
"""

from CM1utils import *

#%% Overview plotting - dbz and thpert

fp = 'C:/Users/mschne28/Documents/cm1out/warmbub_p3_250m/'
# fp = 'C:/Users/mschne28/Documents/cm1r21.1/run/'
fn = np.linspace(1,29,8)


titlestr = "P3 triple-moment, dx=250m"


figsave = False

fig,ax = plt.subplots(2, 4, figsize=(9.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
fig1,ax1 = plt.subplots(2, 4, figsize=(9.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')

for f in fn:
    ds = nc.Dataset(fp+f"cm1out_{f:06.0f}.nc")
    time = ds.variables['time'][:].data[0]
    xh = ds.variables['xh'][:].data
    yh = ds.variables['yh'][:].data
    zh = ds.variables['zh'][:].data
    iz1 = np.where(zh >= 1)[0][1]
    iz2 = np.where(zh >= 2)[0][1]
    iz3 = np.where(zh >= 3)[0][1]
    
    
    dbz = ds.variables['dbz'][:].data[0,0,:,:]
    # cref = ds.variables['cref'][:].data[0,:,:]
    thpert = ds.variables['th'][:].data[0,0,:,:] - ds.variables['th0'][:].data[0,0,:,:]
    uinterp = ds.variables['uinterp'][:].data[0,0:iz1,:,:]
    vinterp = ds.variables['vinterp'][:].data[0,0:iz1,:,:]
    winterp = ds.variables['winterp'][:].data[0,0:iz2,:,:]
    zvort = ds.variables['zvort'][:].data[0,0:iz1,:,:]
    
    ### NSSL scheme
    # thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] - 
    #             (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] + 
    #              ds.variables['qi'][:].data[0,0,:,:] + ds.variables['qs'][:].data[0,0,:,:] + 
    #              ds.variables['qg'][:].data[0,0,:,:] + ds.variables['qhl'][:].data[0,0,:,:]))
    ### P3 scheme
    thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] - 
                (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] + 
                 ds.variables['qi'][:].data[0,0,:,:]))
    thr0 = ds.variables['th0'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,:,:])
    thrpert = thr - thr0
    del thr,thr0
    
    u_gr = uinterp + ds.variables['umove'][:].data[0]
    v_gr = vinterp + ds.variables['vmove'][:].data[0]
    
    ds.close()
    
    
    xl = [-150,150]
    yl = [-150,150]
    
    # xl = [-100,50]
    # yl = [-75,75]
    
    
    n = (f-fn[0])/(fn[1]-fn[0])
    i = int(np.floor(n/4))
    j = int(np.mod(n,4))
    
    if j == 3:
        cb_flag = True
    else:
        cb_flag = False
    
    
    plot_contourf(xh, yh, np.ma.masked_array(dbz, dbz<0.1), 'dbz', ax[i,j], levels=np.linspace(0,70,15),
                  datalims=[0,70], xlims=xl, ylims=yl, cmap='HomeyerRainbow', cbar=cb_flag)
    ax[i,j].contour(xh, yh, np.max(winterp, axis=0), levels=[5,10], colors='k', linestyles='-', linewidths=1)
    # ax[i,j].contour(xh, yh, np.min(winterp, axis=0), levels=[-5], colors='b', linestyles='-', linewidths=1)
    ax[i,j].set_title(f"t = {time:.0f} s")
    fig.suptitle(f"Sfc dbz + max 0-2 km w ({titlestr})")
    if (n==len(fn)-1) & (figsave):
        fig.savefig(fp+f"figs/dbz.png", dpi=300)
    
    
    qix = 60
    
    plot_contourf(xh, yh, thpert, 'thpert', ax1[i,j], levels=np.linspace(-15,15,31),
                  datalims=[-15,15], xlims=xl, ylims=yl, cmap='balance', cbar=cb_flag)
    ax1[i,j].contour(xh, yh, np.max(zvort, axis=0), levels=[0.02], colors='k', linestyles='-', linewidths=1)
    ax1[i,j].quiver(xh[::qix], yh[::qix], u_gr[0,::qix,::qix], v_gr[0,::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
    ax1[i,j].set_title(f"t = {time:.0f} s")
    fig1.suptitle(f"Sfc thpert + sfc wind + max 0-1 km zeta=0.02 s$^{{-1}}$ ({titlestr})")
    if (n==len(fn)-1) & (figsave):
        fig1.savefig(fp+f"figs/thpert.png", dpi=300)





















