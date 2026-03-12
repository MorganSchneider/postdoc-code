# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 09:57:54 2025

@author: mschne28
"""

from CM1utils import *
import matpy.calc as mc
from metpy.units import units

#%% Overview plotting - dbz and thpert

fp = 'C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_250m/'

fn = np.linspace(1,37,10)
ncols = 5


if 'semislip' in fp:
    bbc = 'Semi-slip'
    sim = 'SEMISLIP'
elif 'freeslip' in fp:
    bbc = 'Free-slip'
    sim = 'FREESLIP'
elif 'noslip' in fp:
    bbc = 'No-slip'
    sim = 'NOSLIP'

titlestr = f"{bbc}, P3 3mom, dx=250m"
# titlestr = "New P3 -- Fir, modded"


figsave = False


# fig,ax = plt.subplots(2, ncols, figsize=(9.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
# fig1,ax1 = plt.subplots(2, ncols, figsize=(9.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
# fig2,ax2 = plt.subplots(2, ncols, figsize=(10,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')

# fig,ax = plt.subplots(2, ncols, figsize=(11.75,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
fig1,ax1 = plt.subplots(2, ncols, figsize=(12,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
fig2,ax2 = plt.subplots(2, ncols, figsize=(12.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')

for f in fn:
    ds = nc.Dataset(fp+f"cm1out_{f:06.0f}.nc")
    time = ds.variables['time'][:].data[0]
    xh = ds.variables['xh'][:].data
    yh = ds.variables['yh'][:].data
    zh = ds.variables['zh'][:].data
    iz1 = np.where(zh>1)[0][0]
    iz2 = np.where(zh>2)[0][0]
    iz3 = np.where(zh>3)[0][0]
    
    
    # dbz = ds.variables['dbz'][:].data[0,0,:,:]
    # winterp = ds.variables['winterp'][:].data[0,0:iz2,:,:]
    zvort = ds.variables['zvort'][:].data[0,0:iz1,:,:]
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
    # thrpert = thr - thr0
    # del thr,thr0
    del2 = mc.laplacian(thrpert*units.K, deltas=(250*units.m, 250*units.m))
    del2thp = del2.magnitude
    del2units = del2.units
    
    ds.close()
    
    
    xl = [-150,150]
    yl = [-150,150]
    
    # xl = [-100,100]
    # yl = [-100,100]
    
    
    n = (f-fn[0])/(fn[1]-fn[0])
    i = int(np.floor(n/ncols))
    j = int(np.mod(n,ncols))
    
    if j == ncols-1:
        cb_flag = True
    else:
        cb_flag = False
    
    
    qix = 60
    
    
    # plot_contourf(xh, yh, np.ma.masked_array(dbz, dbz<0.1), 'dbz', ax[i,j], levels=np.linspace(0,70,15),
    #               datalims=[0,70], xlims=xl, ylims=yl, cmap='HomeyerRainbow', cbar=cb_flag)
    # ax[i,j].contour(xh, yh, np.max(winterp, axis=0), levels=[5,10], colors=['dimgray','k'], linestyles='-', linewidths=[0.75,0.75])
    # if n == 0:
    #     l1, = ax[0,0].plot([190,200], [190,200], color='dimgray', linewidth=0.75)
    #     l2, = ax[0,0].plot([190,200], [190,200], '-k', linewidth=0.75)
    #     ax[0,0].legend(handles=[l1,l2], labels=['w=5 m/s','w=10 m/s'], loc='upper right', fontsize=10)
    # ax[i,j].set_title(f"t = {time:.0f} s")
    # # fig.suptitle(f"Sfc dbz + max 0-2 km w ({titlestr})")
    # fig.suptitle(f"{sim}")
    # if (n==len(fn)-1) & (figsave):
    #     fig.savefig(fp+f"figs/dbz.png", dpi=300)
    
    
    
    plot_contourf(xh, yh, thrpert, 'thpert', ax1[i,j], levels=np.linspace(-12,12,25),
                  datalims=[-12,12], xlims=xl, ylims=yl, cmap='balance', cbar=cb_flag)
    ax1[i,j].contour(xh, yh, np.max(zvort, axis=0), levels=[0.03], colors='r', linestyles='-', linewidths=1)
    ax1[i,j].quiver(xh[::qix], yh[::qix], u_gr[::qix,::qix], v_gr[::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
    if n == 0:
        l3, = ax1[0,0].plot([190,200], [190,200], '-r', linewidth=1)
        ax1[0,0].legend(handles=[l3], labels=["\u03B6=0.03 s$^{-1}$"], loc='upper right', fontsize=10)
    ax1[i,j].set_title(f"t = {time:.0f} s")
    # fig1.suptitle(f"Sfc thrpert + sfc wind + max 0-1 km zeta=0.025 s$^{{-1}}$ ({titlestr})")
    fig1.suptitle(f"{sim}")
    if (n==len(fn)-1) & (figsave):
        fig1.savefig(fp+f"figs/thrpert.png", dpi=300)
    
    
    
    c = ax2[i,j].contourf(xh, yh, del2thp, levels=np.linspace(-5e-5,5e-5,21), vmin=-5e-5, vmax=5e-5, cmap='balance', antialiased=True)
    c.set_edgecolor('face')
    if cb_flag:
        cb = plt.colorbar(c, ax=ax2[i,j], extend='both')
        cb.set_label("\u25BD$^2$\u03B8' (K m$^{-2}$)", fontsize=11)
        cb.formatter.set_powerlimits((0,0))
        cb.set_ticks(np.linspace(-5e-5,5e-5,11))
    ax2[i,j].set_xlim(xl)
    ax2[i,j].set_ylim(yl)
    ax2[i,j].quiver(xh[::qix], yh[::qix], u_gr[::qix,::qix], v_gr[::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
    ax2[i,j].set_title(f"t = {time:.0f} s")
    fig2.suptitle(f"{sim} -- \u25BD$^2$\u03B8'")
    if (n==len(fn)-1) & (figsave):
        fig2.savefig(fp+f"figs/thp_laplacian.png", dpi=300)
    
    





#%% Plot swaths

fp = 'C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_250m/'

if 'semislip' in fp:
    bbc = 'Semi-slip'
    sim = 'SEMISLIP'
elif 'freeslip' in fp:
    bbc = 'Free-slip'
    sim = 'FREESLIP'
elif 'noslip' in fp:
    bbc = 'No-slip'
    sim = 'NOSLIP'

# dx = fp[47:-2]

# if 'wk' in fp:
#     bs = 'WK profile'

titlestr = f"{bbc}, WK profile, dx=250"

ds = nc.Dataset(fp+"cm1out_000009.nc")
time = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
# dbz = ds.variables['dbz'][:].data[0,0,:,:]
# wspd = np.sqrt(ds.variables['uinterp'][:].data[0,0,:,:]**2 + ds.variables['vinterp'][:].data[0,0,:,:]**2)

# hail = ds.variables['hail'][:].data[0,:,:]
# hail2 = ds.variables['hail2'][:].data[0,:,:]
sws = ds.variables['sws'][:].data[0,:,:] #max sfc wind speed
svs = ds.variables['svs'][:].data[0,:,:] #max sfc zeta
sps = ds.variables['sps'][:].data[0,:,:] #min sfc p'
sus = ds.variables['sus'][:].data[0,:,:] #max 5km w
shs = ds.variables['shs'][:].data[0,:,:] #max integrated UH
sus2 = ds.variables['sus2'][:].data[0,:,:]
svs2 = ds.variables['svs2'][:].data[0,:,:]
ds.close()


xl = [-150,150]
yl = [-150,150]


figsave = False

# print(np.max(wspd))



                      

# fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
# c = plot_contourf(xh, yh, sws, 'wspd', ax,
#               levels=np.linspace(0,40,21), datalims=[0,40], 
#               xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max')
# # ax.contour(xh, yh, svs, levels=[0.01,0.02], colors=[''])
# ax.contour(xh, yh, sws, levels=[26,33], colors='k', linewidths=[0.5,1.25])
# # ax.set_title(f"7-h sfc wind swath ({titlestr})", fontsize=16)
# ax.set_title(f"{sim} - Surface wind swath", fontsize=16)
# l1, = ax.plot([190,200], [190,200], '-k', linewidth=0.5)
# l2, = ax.plot([190,200], [190,200], '-k', linewidth=1)
# ax.legend(handles=[l1,l2], labels=['Severe (26 m/s)','Sig. Severe (33 m/s)'], loc='upper right', fontsize=12)
# if figsave:
#     plt.savefig(fp+"figs/wspd_swath.png", dpi=300)


# fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
# plot_contourf(xh, yh, np.ma.masked_array(svs, svs<0.001), 'zvort', ax,
#               levels=np.linspace(0,0.05,21), datalims=[0,0.05],
#               xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,0.05,11), extend='max')
# # ax.set_title(f"7-h sfc zeta swath ({titlestr})", fontsize=16)
# ax.set_title(f"{sim} - Surface vorticity swath", fontsize=16)
# # ax.contour(xh, yh, sws, levels=[26,33], colors='k', linewidths=[0.5,1])
# # l1, = ax.plot([0,0], [190,200], '-k', linewidth=0.5)
# # l2, = ax.plot([0,0], [190,200], '-k', linewidth=1)
# # ax.legend(handles=[l1,l2], labels=['Severe (26 m/s)','Sig. Severe (33 m/s)'], loc='upper right', fontsize=12)
# ax.contour(xh, yh, sus, levels=[20,40], colors='k', linewidths=[0.75,1.5])
# l1, = ax.plot([190,200], [190,200], '-k', linewidth=0.75)
# l2, = ax.plot([190,200], [190,200], '-k', linewidth=1.5)
# ax.legend(handles=[l1,l2], labels=['w=20 m/s','w=40 m/s'], loc='upper right', fontsize=12)
# if figsave:
#     plt.savefig(fp+"figs/zeta_swath.png", dpi=300)


fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, np.ma.masked_array(sus, sus<1), 'w', ax,
              levels=np.linspace(0,40,21), datalims=[0,40],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max')
ax.contour(xh, yh, svs, levels=[0.04], colors='k', linewidths=[0.75])
l1, = ax.plot([190,200], [190,200], '-k', linewidth=0.75)
l2, = ax.plot([190,200], [190,200], '-k', linewidth=1)
ax.legend(handles=[l1], labels=["\u03B6=0.04 s$^{-1}$"], loc='upper right', fontsize=12)
ax.set_title(f"{sim} - 5-km updraft swath", fontsize=16)
if figsave:
    plt.savefig(fp+"figs/updraft_swath.png", dpi=300)


fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, np.ma.masked_array(sus2, sus2<1), 'w', ax,
              levels=np.linspace(0,40,21), datalims=[0,40],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max')
ax.contour(xh, yh, svs2, levels=[0.04], colors='k', linewidths=[0.75])
l1, = ax.plot([190,200], [190,200], '-k', linewidth=0.75)
l2, = ax.plot([190,200], [190,200], '-k', linewidth=1)
ax.legend(handles=[l1], labels=["\u03B6=0.04 s$^{-1}$"], loc='upper right', fontsize=12)
ax.set_title(f"{sim} - Translated 5-km updraft swath", fontsize=16)


# fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
# plot_contourf(xh, yh, np.ma.masked_array(shs, shs<10), 'uh', ax,
#               levels=np.linspace(0,2000,21), datalims=[0,2000],
#               xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,2000,11), extend='max')
# ax.set_title(f"7-h UH swath ({titlestr})", fontsize=16)
# if figsave:
#     plt.savefig(fp+"figs/helicity_swath.png", dpi=300)


# fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
# plot_contourf(xh, yh, np.ma.masked_array(sps/100, sps/100>-0.1), 'prspert', ax,
#               levels=np.linspace(-5,0,21), datalims=[-5,0],
#               xlims=xl, ylims=yl, cmap='Blues_r', cbfs=14, cbticks=np.linspace(-5,0,11), extend='min')
# ax.set_title(f"7-h sfc p' swath ({titlestr})", fontsize=16)
# if figsave:
#     plt.savefig(fp+"figs/prspert_swath.png", dpi=300)


plt.show()




# Hourly statistics?






#%% Multi-panel swaths for CWE?

# Freeslip
ds = nc.Dataset('C:/Users/mschne28/Documents/cm1out/cwe/freeslip_wk_250m/cm1out_000029.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
sws_fs = ds.variables['sws'][:].data[0,:,:] #max sfc wind
svs_fs = ds.variables['svs'][:].data[0,:,:] #max sfc vort
sus_fs = ds.variables['sus'][:].data[0,:,:] #max 5km updraft
shs_fs = ds.variables['shs'][:].data[0,:,:] #max integrated UH
ds.close()


# Semislip
ds = nc.Dataset('C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_250m/cm1out_000029.nc')
sws_ss = ds.variables['sws'][:].data[0,:,:]
svs_ss = ds.variables['svs'][:].data[0,:,:]
sus_ss = ds.variables['sus'][:].data[0,:,:]
shs_ss = ds.variables['shs'][:].data[0,:,:]
ds.close()

# Noslip
ds = nc.Dataset('C:/Users/mschne28/Documents/cm1out/cwe/noslip_wk_250m/cm1out_000029.nc')
sws_ns = ds.variables['sws'][:].data[0,:,:]
svs_ns = ds.variables['svs'][:].data[0,:,:]
sus_ns = ds.variables['sus'][:].data[0,:,:]
shs_ns = ds.variables['shs'][:].data[0,:,:]
ds.close()


xl = [-150,150]
yl = [-150,150]


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


plt.show()

#%% plot cm1 stats time series

fp = 'C:/Users/mschne28/Documents/cm1out/cwe/noslip_wk_250m/'


if 'semislip' in fp:
    bbc = 'Semi-slip'
    sim = 'SEMISLIP'
elif 'freeslip' in fp:
    bbc = 'Free-slip'
    sim = 'FREESLIP'
elif 'noslip' in fp:
    bbc = 'No-slip'
    sim = 'NOSLIP'

# dx = fp[47:-2]

titlestr = f"{bbc}, WK profile, dx=250m"

# ds = nc.Dataset(fp+'cm1out_000001.nc')
# prs0 = ds.variables['prs0'][:].data[0,0,0,0]
# ds.close()

ds = nc.Dataset(fp+f"cm1out_stats.nc")
time = ds.variables['mtime'][:].data
wmax500 = ds.variables['wmax500'][:].data #max w at 500 m
wmin500 = ds.variables['wmin500'][:].data #min w at 500 m
wmax1000 = ds.variables['wmax1000'][:].data #max w at 1000 m
wmin1000 = ds.variables['wmin1000'][:].data #min w at 1000 m
wmax2500 = ds.variables['wmax2500'][:].data #max w at 2500 m
wmin2500 = ds.variables['wmin2500'][:].data #min w at 2500 m
wmax5000 = ds.variables['wmax5000'][:].data #max w at 5000 m
wmin5000 = ds.variables['wmin5000'][:].data #min w at 5000 m
swspmax = ds.variables['swspmax'][:].data #max sfc wspd
divmax = ds.variables['divmax'][:].data #max 3d divergence
divmin = ds.variables['divmin'][:].data #max 3d convergence
vortsfc = ds.variables['vortsfc'][:].data #max sfc vort
vort1km = ds.variables['vort1km'][:].data #max 1km vort
vort2km = ds.variables['vort2km'][:].data #max 2km vort
vort3km = ds.variables['vort3km'][:].data #max 3km vort
vort4km = ds.variables['vort4km'][:].data #max 4km vort
vort5km = ds.variables['vort5km'][:].data #max 5km vort
tkemax = ds.variables['tkemax'][:].data #max subgrid tke
# sprsmax = ds.variables['sprsmax'][:].data - prs0 #max sfc p'
# sprsmin = ds.variables['sprsmin'][:].data - prs0 #min sfc p'
ds.close()





# fig,ax = plt.subplots(2, 1, figsize=(10,8), sharex=True, layout='constrained')

# l1,= ax[0].plot(time, movmean(wmax500,5), 'lightcoral', linewidth=2)
# l2,= ax[0].plot(time, movmean(wmax1000,5), 'r', linewidth=2)
# l3,= ax[0].plot(time, movmean(wmax2500,5), 'maroon', linewidth=2)
# l4,= ax[0].plot(time, movmean(wmax5000,5), 'darkviolet', linewidth=2)
# ax[0].set_xlim([0,25200])
# ax[0].set_ylim([0,50]) #24 without 5000m, 50 with 5000m
# ax[0].set_ylabel('max w (m/s)', fontsize=14)
# ax[0].tick_params(axis='both', labelsize=12)
# ax[0].set_title(f"max w ({fp[35:-1]})", fontsize=16)
# ax[0].legend(handles=[l1,l2,l3,l4], labels=['500 m', '1000 m', '2500 m', '5000 m'],
#              loc='upper left', fontsize=14)
# ax[0].grid(visible=True, which='major', color='darkgray', linestyle='-')
# ax[0].grid(visible=True, which='minor', color='lightgray', linestyle='-')
# ax[0].yaxis.set_major_locator(MultipleLocator(10)) #4 without 5000m, 10 with 5000m
# ax[0].yaxis.set_minor_locator(MultipleLocator(5)) #2 without 5000m, 5 with 5000m

# l1,= ax[1].plot(time, movmean(wmin500,5), 'deepskyblue', linewidth=2)
# l2,= ax[1].plot(time, movmean(wmin1000,5), 'b', linewidth=2)
# l3,= ax[1].plot(time, movmean(wmin2500,5), 'navy', linewidth=2)
# l4,= ax[1].plot(time, movmean(wmin5000,5), 'k', linewidth=2)
# ax[1].set_xlim([0,25200])
# ax[1].set_ylim([-20,0])
# ax[1].set_xlabel('Time (s)', fontsize=14)
# ax[1].set_ylabel('min w (m/s)', fontsize=14)
# ax[1].tick_params(axis='both', labelsize=12)
# ax[1].set_title(f"min w ({fp[35:-1]})", fontsize=16)
# ax[1].legend(handles=[l1,l2,l3,l4], labels=['500 m', '1000 m', '2500 m', '5000 m'],
#              loc='lower left', fontsize=14)
# ax[1].grid(visible=True, which='major', color='darkgray', linestyle='-')
# ax[1].grid(visible=True, which='minor', color='lightgray', linestyle='-')
# ax[1].xaxis.set_major_locator(MultipleLocator(3600))
# ax[1].xaxis.set_minor_locator(MultipleLocator(900))
# ax[1].yaxis.set_major_locator(MultipleLocator(4))
# ax[1].yaxis.set_minor_locator(MultipleLocator(2))

# plt.show()


figsave = False



fig,ax = plt.subplots(2, 1, figsize=(10,8), sharex=True, layout='constrained')

# l1,= ax[0].plot([0,25200], [20,20], 'tab:orange', linewidth=1.25)
# l2,= ax[0].plot([0,25200], [26,26], 'red', linewidth=1.25)
# l3,= ax[0].plot([0,25200], [33,33], 'darkviolet', linewidth=1.25)
# ax[0].plot(time, swspmax, 'k', linewidth=2)
# ax[0].set_xlim([0,25200])
# ax[0].set_ylim([0,62])
# # ax[0].set_xlabel('Time (s)', fontsize=14)
# ax[0].set_ylabel('Wind speed (m/s)', fontsize=14)
# ax[0].tick_params(axis='both', labelsize=12)
# # ax[0].set_title(f"Max lml wind speed ({titlestr})", fontsize=16)
# ax[0].set_title(f"{sim} - Max. wind speed", fontsize=16)
# ax[0].grid(visible=True, which='major', color='darkgray', linestyle='-')
# ax[0].grid(visible=True, which='minor', color='lightgray', linestyle='-')
# ax[0].yaxis.set_major_locator(MultipleLocator(10))
# ax[0].yaxis.set_minor_locator(MultipleLocator(5))
# ax[0].legend(handles=[l1,l2,l3], labels=['Sub-severe', 'Severe', 'Sig. Severe'],
#              loc='upper left', fontsize=14)

l1,= ax[0].plot(time, movmean(wmax500,5), 'lightcoral', linewidth=2)
l2,= ax[0].plot(time, movmean(wmax1000,5), 'r', linewidth=2)
l3,= ax[0].plot(time, movmean(wmax2500,5), 'maroon', linewidth=2)
l4,= ax[0].plot(time, movmean(wmax5000,5), 'darkviolet', linewidth=2)
ax[0].set_xlim([0,25200])
ax[0].set_ylim([0,50]) #24 without 5000m, 50 with 5000m
ax[0].set_ylabel('w (m/s)', fontsize=14)
ax[0].tick_params(axis='both', labelsize=12)
ax[0].set_title(f"{sim} - Max. updraft speed", fontsize=16)
ax[0].legend(handles=[l1,l2,l3,l4], labels=['500 m', '1000 m', '2500 m', '5000 m'],
             loc='upper left', fontsize=14)
ax[0].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[0].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[0].yaxis.set_major_locator(MultipleLocator(10)) #4 without 5000m, 10 with 5000m
ax[0].yaxis.set_minor_locator(MultipleLocator(5)) #2 without 5000m, 5 with 5000m


l1,= ax[1].plot(time, movmean(vortsfc,5), 'lightcoral', linewidth=2)
l2,= ax[1].plot(time, movmean(vort1km,5), 'r', linewidth=2)
l3,= ax[1].plot(time, movmean(vort3km,5), 'maroon', linewidth=2)
ax[1].set_xlim([0,25200])
ax[1].set_ylim([0,0.25])
ax[1].set_xlabel('Time (s)', fontsize=14)
ax[1].set_ylabel("Vorticity (1/s)", fontsize=14)
ax[1].tick_params(axis='both', labelsize=12)
# ax[1].set_title(f"Max zeta ({titlestr})", fontsize=16)
ax[1].set_title(f"{sim} - Max. vertical vorticity")
ax[1].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[1].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[1].xaxis.set_major_locator(MultipleLocator(3600))
ax[1].xaxis.set_minor_locator(MultipleLocator(900))
ax[1].yaxis.set_major_locator(MultipleLocator(0.05))
ax[1].yaxis.set_minor_locator(MultipleLocator(0.01))
ax[1].legend(handles=[l1,l2,l3], labels=['10m','1km','3km'],
             loc='upper left', fontsize=14)

if figsave:
    plt.savefig(fp+'figs/w_zeta_ts.png', dpi=300)

plt.show()






# cm1out_stats variables
# time/mtime: time and model time (s since beginning of simulation)

# wmax/wmin/wmax500/1000/2500/5000/10k: max and min w in domain and at 500/1000/2500/5000/10000 m
# zwmax/zwmin: level of max and min w
# umax/umin/vmax/vmin: max and min u and v
# sumax/sumin/svmax/svmin: max and min u and v at LML
# wspmax/wspmin/swsp/wsp10: max and min wind speed in domain and at LML and 10 m
# zwspmax/zwspmin: level of max and min wind speed
# divmax/divmin: max and min 3d divergence
# vortsfc/1km/etc: max zeta at LML/1/2/3/4/5 km

# ppimin/ppimin/ppmax/ppmin: max and min pi and prs pert
# sprsmax/sprsmin/psfcmax/psfcmin: max and min pressure at LML and surface
# thpmax/thpmin/sthpmax/sthpmin: max and min theta pert in domain and at LML
# themax/themin/sthemax/sthemin: max and min theta-e below 10 km and at LML

# maxqv/minqv/qc/qr/qi/qs/qg/qhl: max and min mixing ratios (water vapor, cloud water, rain, ice, snow, graupel, hail)
# maxccn/minccn/ccw/crw/cci/csw/chw/chl: max and min number concentrations (ccn, cloud water, rain, ice, snow, graupel, hail)
# maxvhw/minvhw/maxvhl/minvhl: max and min graupel and hail density
# pratemax/pratemin: max and min surface precip rate
# rhmax/rhmin/rhimax/rhimin: max and min RH wrt liquid and ice
# qctop/qcbot: max cloud top height and min cloud base height

# tkemax/tkemin: max and min subgrid TKE
# hpblmax/hpblmin: max diagosed PBL depth
# kmhmax/kmhmin/kmvmax/kmvmin: max and min horiz and vert eddy viscosity for momentum
# khhmax/khhmin/khvmax/khvmin: max and min horiz and vert eddy diffusivity for scalars
# cflmax/kshmax/ksvmax: max Courant number and horiz/vert K stability factor



#%% Plots for CWE 2026 extended abstract

fp = 'C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_250m/'

# ds = nc.Dataset(fp+'cm1out_000001.nc')
# prs0 = ds.variables['prs0'][:].data[0,0,0,0]
# ds.close()

ds = nc.Dataset(fp+f"cm1out_stats.nc")
time = ds.variables['mtime'][:].data
wmax500 = ds.variables['wmax500'][:].data #max w at 500 m
wmax1000 = ds.variables['wmax1000'][:].data #max w at 1000 m
wmax2500 = ds.variables['wmax2500'][:].data #max w at 2500 m
wmax5000 = ds.variables['wmax5000'][:].data #max w at 5000 m
swspmax = ds.variables['swspmax'][:].data #max lml wspd
vortsfc = ds.variables['vortsfc'][:].data #max lml vort
vort1km = ds.variables['vort1km'][:].data #max 1km vort
vort2km = ds.variables['vort2km'][:].data #max 2km vort
vort3km = ds.variables['vort3km'][:].data #max 3km vort
ds.close()


fp = 'C:/Users/mschne28/Documents/cm1out/cwe/freeslip_wk_250m/'
ds = nc.Dataset(fp+f"cm1out_stats.nc")
wmax500_fs = ds.variables['wmax500'][:].data #max w at 500 m
wmax1000_fs = ds.variables['wmax1000'][:].data #max w at 1000 m
wmax2500_fs = ds.variables['wmax2500'][:].data #max w at 2500 m
wmax5000_fs = ds.variables['wmax5000'][:].data #max w at 5000 m
swspmax_fs = ds.variables['swspmax'][:].data #max sfc wspd
vortsfc_fs = ds.variables['vortsfc'][:].data #max sfc vort
vort1km_fs = ds.variables['vort1km'][:].data #max 1km vort
vort2km_fs = ds.variables['vort2km'][:].data #max 2km vort
vort3km_fs = ds.variables['vort3km'][:].data #max 3km vort
ds.close()


fp = 'C:/Users/mschne28/Documents/cm1out/cwe/noslip_wk_250m/'
ds = nc.Dataset(fp+f"cm1out_stats.nc")
wmax500_ns = ds.variables['wmax500'][:].data #max w at 500 m
wmax1000_ns = ds.variables['wmax1000'][:].data #max w at 1000 m
wmax2500_ns = ds.variables['wmax2500'][:].data #max w at 2500 m
wmax5000_ns = ds.variables['wmax5000'][:].data #max w at 5000 m
swspmax_ns = ds.variables['swspmax'][:].data #max sfc wspd
vortsfc_ns = ds.variables['vortsfc'][:].data #max sfc vort
vort1km_ns = ds.variables['vort1km'][:].data #max 1km vort
vort2km_ns = ds.variables['vort2km'][:].data #max 2km vort
vort3km_ns = ds.variables['vort3km'][:].data #max 3km vort
ds.close()


figsave = True



fig,ax = plt.subplots(3, 1, figsize=(10,9), sharex=True, layout='constrained')

l1,= ax[0].plot(time, movmean(vortsfc_fs,5), 'k', linewidth=2)
l2,= ax[0].plot(time, movmean(vortsfc_ns,5), 'dodgerblue', linewidth=2)
l3,= ax[0].plot(time, movmean(vortsfc,5), 'crimson', linewidth=2)
# l1,= ax[0].plot(time, vortsfc_fs, 'k', linewidth=2)
# l2,= ax[0].plot(time,vortsfc_ns, 'dodgerblue', linewidth=2)
# l3,= ax[0].plot(time, vortsfc, 'crimson', linewidth=2)
ax[0].set_xlim([0,25200])
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
l6,= ax[1].plot(time, movmean(vort1km,5), 'crimson', linewidth=2)
ax[1].set_xlim([0,25200])
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
l9,= ax[2].plot(time, movmean(vort3km,5), 'crimson', linewidth=2)
ax[2].set_xlim([0,25200])
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
    plt.savefig('C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_250m/figs/zeta_all_timeseries.png', dpi=300)
# plt.show()

#%%


fig,ax = plt.subplots(3, 1, figsize=(10,9), sharex=True, layout='constrained')

l1,= ax[0].plot(time, movmean(wmax1000_fs,5), 'k', linewidth=2)
l2,= ax[0].plot(time, movmean(wmax1000_ns,5), 'dodgerblue', linewidth=2)
l3,= ax[0].plot(time, movmean(wmax1000,5), 'crimson', linewidth=2)
ax[0].set_xlim([0,25200])
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
l6,= ax[1].plot(time, movmean(wmax2500,5), 'crimson', linewidth=2)
ax[1].set_xlim([0,25200])
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
l9,= ax[2].plot(time, movmean(wmax5000,5), 'crimson', linewidth=2)
ax[2].set_xlim([0,25200])
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
    plt.savefig('C:/Users/mschne28/Documents/cm1out/cwe/semislip_wk_250m/figs/w_all_timeseries.png', dpi=300)

# plt.show()



#%% Zooming into the TLV in freeslip? (OLD P3 SCHEME)

fp = 'C:/Users/mschne28/Documents/cm1out/freeslip_wk_250m/'
# fp = 'C:/Users/mschne28/Documents/cm1r21.1/run/'
fn = 25


if 'semislip' in fp:
    bbc = 'Semi-slip'
elif 'freeslip' in fp:
    bbc = 'Free-slip'
elif 'noslip' in fp:
    bbc = 'No-slip'

titlestr = f"{bbc}, WK profile"







xl = [-60,-10]
yl = [-110,-60]


ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
time = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['zh'][:].data
iz1 = np.where(zh>1)[0][0]
iz2 = np.where(zh>2)[0][0]
iz3 = np.where(zh>3)[0][0]
iz5 = np.where(zh>5)[0][0]
ix1 = np.argmin(np.abs(xh-xl[0]))
ix2 = np.argmin(np.abs(xh-xl[1]))
iy1 = np.argmin(np.abs(yh-yl[0]))
iy2 = np.argmin(np.abs(yh-yl[1]))
ix = slice(ix1,ix2)
iy = slice(iy1,iy2)



xh = xh[ix]
yh = yh[iy]
zh = zh[0:iz5]

dbz = ds.variables['dbz'][:].data[0,0,iy,ix]
thpert = ds.variables['th'][:].data[0,0,iy,ix] - ds.variables['th0'][:].data[0,0,iy,ix]
uinterp = ds.variables['uinterp'][:].data[0,0:iz1,iy,ix]
vinterp = ds.variables['vinterp'][:].data[0,0:iz1,iy,ix]
winterp = ds.variables['winterp'][:].data[0,0:iz2,iy,ix]
zvort = ds.variables['zvort'][:].data[0,0:iz1,iy,ix]

ic = np.where(np.max(zvort, axis=0) == np.max(np.max(zvort, axis=0)))
ixc = slice(ic[1][0]-40, ic[1][0]+41)
iyc = ic[0][0]

zvortc = ds.variables['zvort'][:].data[0,0:iz5,iy,ix][:,iyc,ixc]
winterpc = ds.variables['winterp'][:].data[0,0:iz5,iy,ix][:,iyc,ixc]
thpertc = ds.variables['th'][:].data[0,0:iz5,iy,ix][:,iyc,ixc] - ds.variables['th0'][:].data[0,0:iz5,iy,ix][:,iyc,ixc]

### P3 scheme
thr = ds.variables['th'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,0,iy,ix] - 
            (ds.variables['qc'][:].data[0,0,iy,ix] + ds.variables['qr'][:].data[0,0,iy,ix] + 
             ds.variables['qi'][:].data[0,0,iy,ix]))
# thr = ds.variables['th'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,0,iy,ix] - #for updated P3 scheme
#             (ds.variables['qc'][:].data[0,0,iy,ix] + ds.variables['qr'][:].data[0,0,iy,ix] + 
#              ds.variables['qi1'][:].data[0,0,iy,ix] + ds.variables['qi2'][:].data[0,0,iy,ix] +
#              ds.variables['qi3'][:].data[0,0,iy,ix] + ds.variables['qi4'][:].data[0,0,iy,ix]))
thr0 = ds.variables['th0'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,iy,ix])
thrpert = thr - thr0
del thr,thr0

u_gr = uinterp + ds.variables['umove'][:].data[0]
v_gr = vinterp + ds.variables['vmove'][:].data[0]

ds.close()


#%%

figsave = False



xleg = [200,201]
yleg = [200,201]

# zvlevs = [0.025, 0.05, 0.1]
# zvlabs = [f"\u03B6={zvlevs[i]} s$^{-1}$" for i in range(len(zvlevs))]
# zvcols = ['k', 'k', 'k']
# zvlws = [1, 1.5, 2]
# wlevs = [5, 10]
# wlabs = [f"w={wlevs[i]} m/s" for i in range(len(wlevs))]
# wcols = ['xkcd:slate', 'xkcd:slate']
# wlws = [0.5, 1]
# labs = list(np.append(zvlabs, wlabs))
# cols = list(np.append(zvcols, wcols))
# lws = list(np.append(zvlws, wlws))
# l = []


fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, np.ma.masked_array(dbz, dbz<0.1), 'dbz', ax, levels=np.linspace(0,70,15),
              datalims=[0,70], xlims=xl, ylims=yl, cmap='HomeyerRainbow')
ax.contour(xh, yh, np.max(winterp, axis=0), levels=[5,10], colors='xkcd:slate', linestyles='-', linewidths=[0.5,1])
ax.contour(xh, yh, np.max(zvort, axis=0), levels=[0.025, 0.05, 0.1], colors='k', linestyles='-', linewidths=[1,1.5,2])
l1, = ax.plot(xleg, yleg, 'k', linewidth=1, label="\u03B6=0.025 s$^{-1}$")
l2, = ax.plot(xleg, yleg, 'k', linewidth=1.5, label="\u03B6=0.05 s$^{-1}$")
l3, = ax.plot(xleg, yleg, 'k', linewidth=2, label="\u03B6=0.1 s$^{-1}$")
l4, = ax.plot(xleg, yleg, c='xkcd:slate', linewidth=0.5, label="w=5 m/s")
l5, = ax.plot(xleg, yleg, c='xkcd:slate', linewidth=1, label="w=10 m/s")
ax.legend(handles=[l1,l2,l3,l4,l5], ncol=2, loc='upper left', fontsize=10)
ax.set_title(f"Sfc dbz, max 0-2 km w, max 0-1 km \u03B6, t={time:.0f} s\n ({titlestr}) ")
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
if figsave:
    fig.savefig(fp+f"figs/dbz.png", dpi=300)


qix = 8


# zvlevs = [0.025, 0.05, 0.1]
# zvlabs = [f"\u03B6={zvlevs[i]} s$^{-1}$" for i in range(len(zvlevs))]
# zvcols = ['k','k','r']
# zvlws = [1,2,2]
# l = []


fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, thrpert, 'thpert', ax, levels=np.linspace(-12,12,25),
              datalims=[-12,12], xlims=xl, ylims=yl, cmap='balance')
ax.contour(xh, yh, np.max(zvort, axis=0), levels=[0.025,0.05,0.1], colors=['k','k','r'], linestyles='-', linewidths=[1,2,2])
ax.quiver(xh[::qix], yh[::qix], u_gr[0,::qix,::qix], v_gr[0,::qix,::qix], color='k', scale=500, width=0.003, pivot='middle')
l1, = ax.plot(xleg, yleg, 'k', linewidth=1, label="\u03B6=0.025 s$^{-1}$")
l2, = ax.plot(xleg, yleg, 'k', linewidth=2, label="\u03B6=0.05 s$^{-1}$")
l3, = ax.plot(xleg, yleg, 'r', linewidth=2, label="\u03B6=0.1 s$^{-1}$")
ax.legend(handles=[l1,l2,l3], loc='upper left', fontsize=10)
# for c,lw in zip(zvcols,zvlws):
#     l1, = ax.plot(xleg, yleg, color=c, linewidth=lw)
#     l = list(np.append(l, l1))
# ax.legend(handles=l, labels=zvlabs, loc='upper left', fontsize=10)
ax.set_title(f"Sfc \u03B8'\u1D68, max 0-1 km \u03B6, t={time:.0f} s\n ({titlestr}) ")
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
if figsave:
    fig.savefig(fp+f"figs/thrpert.png", dpi=300)






#%%

zvlevs = [-0.005, 0.01, 0.025, 0.05, 0.1]
zvlabs = [f"\u03B6={zvlevs[i]} s$^{{-1}}$" for i in range(len(zvlevs))]
zvlws = [0.5, 0.5, 0.75, 1, 2]
zvlss = ['--', '-', '-', '-', '-']
l = []

xlc = [xh[ixc][0], xh[ixc][-1]]


fig,ax = plt.subplots(1, 1, figsize=(9,4), subplot_kw=dict(box_aspect=0.5), layout='constrained')
plot_contourf(xh[ixc], zh, winterpc, 'w', ax, levels=np.linspace(-36,36,37),
              datalims=[-35,35], xlims=xlc, ylims=[0,5], cmap='balance')
ax.contour(xh[ixc], zh, zvortc, levels=zvlevs, colors='k', linestyles=zvlss, linewidths=zvlws)
# l1, = ax.plot(xleg, yleg, 'k', linewidth=0.5, label="\u03B6=0.01 s$^{-1}$")
# l2, = ax.plot(xleg, yleg, 'k', linewidth=0.75, label="\u03B6=0.025 s$^{-1}$")
# l3, = ax.plot(xleg, yleg, 'k', linewidth=1, label="\u03B6=0.05 s$^{-1}$")
# l4, = ax.plot(xleg, yleg, 'k', linewidth=2, label="\u03B6=0.1 s$^{-1}$")
# ax.legend(handles=[l1,l2,l3], loc='upper left', fontsize=10)
for ls,lw in zip(zvlss,zvlws):
    l1, = ax.plot(xleg, yleg, color='k', linestyle=ls, linewidth=lw)
    l = list(np.append(l, l1))
ax.legend(handles=l, labels=zvlabs, loc='upper left', fontsize=8)
ax.set_xlabel('x (km)')
ax.set_ylabel('z (km)')
ax.set_title(f"Cross section through \u03B6-max at t={time:.0f} s\n w shaded, \u03B6 contour ({titlestr}) ")
plt.show()



# wlevs = [-10, -4, 4, 10, 15, 20]
# wlabs = [f"w={wlevs[i]} m/s" for i in range(len(wlevs))]
# wcols = ['b', 'b', 'k', 'k', 'k', 'k']
# wlss = ['--', '--', '-', '-', '-', '-']
# wlws = [1, 0.5, 0.5, 1, 1.5, 2]
# l = []

# fig,ax = plt.subplots(1, 1, figsize=(9,4), subplot_kw=dict(box_aspect=0.5), layout='constrained')
# plot_contourf(xh[ixc], zh, zvortc, 'zvort', ax, levels=np.linspace(-0.15,0.15,31),
#               datalims=[-0.15,0.15], xlims=xlc, ylims=[0,5], cmap='balance')
# # ax.contour(xh[ixc], zh, winterpc, levels=[-10,-4], colors='b', linewidths=[1,0.5])
# # ax.contour(xh[ixc], zh, winterpc, levels=[4,10,15,20], colors='k', linewidths=[0.5,1,1.5,2])
# ax.contour(xh[ixc], zh, winterpc, levels=wlevs, colors=wcols, linestyles=wlss, linewidths=wlws)
# for c,ls,lw in zip(wcols,wlss,wlws):
#     l1, = ax.plot(xleg, yleg, color=c, linestyle=ls, linewidth=lw)
#     l = list(np.append(l,l1))
# ax.legend(handles=l, labels=wlabs, loc='upper left', fontsize=8)
# ax.set_xlabel('x (km)')
# ax.set_ylabel('z (km)')
# ax.set_title(f"Cross section through \u03B6-max at t={time:.0f} s\n \u03B6 shaded, w contour ({titlestr}) ")
# plt.show()


#%%

wlevs = [-6, -3, 3, 6, 10, 15, 20]
wcols = ['b', 'b', 'r', 'r', 'r', 'r', 'r']
wlss = ['--', '--', '--', '--', '-', '-', '-']
wlws = [0.75, 0.5, 0.5, 0.75, 0.75, 1.25, 1.75]
wlabs = [f"w={wlevs[i]} m/s" for i in range(len(wlevs))]

zvlevs = [-0.005, 0.01, 0.05, 0.1]
zvcols = ['k', 'k', 'k', 'k']
zvlws = [0.5, 0.75, 1.25, 2]
zvlss = ['--', '-', '-', '-']
zvlabs = [f"\u03B6={zvlevs[i]} s$^{{-1}}$" for i in range(len(zvlevs))]

labs = list(np.append(zvlabs, wlabs))
l = []

fig,ax = plt.subplots(1, 1, figsize=(9,4), subplot_kw=dict(box_aspect=0.5), layout='constrained')
plot_contourf(xh[ixc], zh, thpertc, 'thpert', ax, levels=np.linspace(-12,12,21),
              datalims=[-12,12], xlims=xlc, ylims=[0,5], cmap=cmocean.cm.balance)

cw = ax.contour(xh[ixc], zh, winterpc, levels=wlevs, colors=wcols, linestyles=wlss, linewidths=wlws)
ax.clabel(cw, inline=True, fmt="%.0f", fontsize=7, inline_spacing=2)
for lab in cw.labelTexts:
    lab.set_rotation(0)

# cz = ax.contour(xh[ixc], zh, zvortc, levels=zvlevs, colors=zvcols, linestyles=zvlss, linewidths=zvlws)
# for c,ls,lw in zip(zvcols,zvlss,zvlws): #zvort contours
#     l1, = ax.plot(xleg, yleg, color=c, linestyle=ls, linewidth=lw)
#     l = list(np.append(l, l1))
# ax.legend(handles=l, labels=zvlabs, loc='upper left', fontsize=8)

ax.set_xlabel('x (km)')
ax.set_ylabel('z (km)')
ax.set_title(f"Cross section through \u03B6-max at t={time:.0f} s\n \u03B8' shaded, w and \u03B6 contour ({titlestr}) ")

plt.show()



#%% Base state soundings

from metpy.plots import SkewT, Hodograph
import metpy.calc as mc
from metpy.units import units
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Effective Shear Algorithm
def effective_layer(p, T, Td, h, height_layer=False):
    '''This function determines the effective inflow layer for a convective sounding.
    
    Input:
      - p: sounding pressure with units
      - T: sounding temperature with units
      - Td: sounding dewpoint temperature with units
      - h: sounding heights with units
      
    Returns:
      - pbot/hbot, ptop/htop: pressure/height of the bottom level, pressure/height of the top level
    '''
    
    pbot = None
    
    for i in range(p.shape[0]):
        prof = mc.parcel_profile(p[i:], T[i], Td[i])
        sbcape, sbcin = mc.cape_cin(p[i:], T[i:], Td[i:], prof)
        if sbcape >= 100 * units('J/kg') and sbcin > -250 * units('J/kg'):
            pbot = p[i]
            hbot = h[i]
            bot_idx = i
            break
    if not pbot:
        return None, None
    
    for i in range(bot_idx+1, p.shape[0]):
        prof = mc.parcel_profile(p[i:], T[i], Td[i])
        sbcape, sbcin = mc.cape_cin(p[i:], T[i:], Td[i:], prof)
        if sbcape < 100 * units('J/kg') or sbcin < -250 * units('J/kg'):
            ptop = p[i]
            htop = h[i]
            break
            
    if height_layer:
        return hbot, htop
    else:
        return pbot, ptop

effl = effective_layer(prs0/100*units.hPa, T0*units.K, Td0*units.K, zh*1000*units.m, height_layer=True)
ebot = effl[0].magnitude
etop = effl[1].magnitude


ds = nc.Dataset('C:/Users/mschne28/Documents/cm1out/cwe/semislip_nssl_500m/cm1out_000001.nc')
zh = ds.variables['zh'][:].data # km
umove = ds.variables['umove'][:].data
vmove = ds.variables['vmove'][:].data
u0 = ds.variables['u0'][:].data[0,:,0,0] + umove # m/s
v0 = ds.variables['v0'][:].data[0,:,0,0] + vmove # m/s
th0 = ds.variables['th0'][:].data[0,:,0,0] # K
qv0 = ds.variables['qv0'][:].data[0,:,0,0] # kg/kg
prs0 = ds.variables['prs0'][:].data[0,:,0,0] # Pa
ds.close()

T0 = th0 * (prs0/100000.)**0.286
e0 = (qv0 * prs0/100) / (0.622+qv0)
Td0 = 243.5 / ((17.67/(np.log(e0/6.112)))-1) + 273.15


# Calculate sounding parameters
bwnd = mc.bunkers_storm_motion(prs0/100*units.hPa, u0*units('m/s'), v0*units('m/s'), zh*1000*units.m)
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

T_parcel = mc.parcel_profile(prs0/100*units.hPa, T0[0]*units.K, Td0[0]*units.K)
# CC = mc.cape_cin(prs0/100*units.hPa, T0*units.K, Td0*units.K, T_parcel)
# cape = CC[0].magnitude
# cin = CC[1].magnitude
MU = mc.most_unstable_cape_cin(prs0/100*units.hPa, T0*units.K, Td0*units.K, height=zh*1000*units.m)
mucape = MU[0].magnitude
mucin = MU[1].magnitude
SB = mc.surface_based_cape_cin(prs0/100*units.hPa, T0*units.K, Td0*units.K)
sbcape = SB[0].magnitude
sbcin = SB[1].magnitude
ML = mc.mixed_layer_cape_cin(prs0/100*units.hPa, T0*units.K, Td0*units.K, height=zh*1000*units.m)
mlcape = ML[0].magnitude
mlcin = ML[1].magnitude
DC = mc.downdraft_cape(prs0/100*units.hPa, T0*units.K, Td0*units.K)
dcape = DC[0].magnitude

lcl = mc.lcl(prs0/100*units.hPa, T0[0]*units.K, Td0[0]*units.K)
plcl = lcl[0].magnitude[0]
ilcl = np.where(prs0/100 <= plcl)[0][0]
zlcl = zh[ilcl]*1000

el = mc.el(prs0/100*units.hPa, T0*units.K, Td0*units.K)
pel = el[0].magnitude
iel = np.where(prs0/100 <= pel)[0][0]
zel = zh[iel]*1000


shr6km = mc.bulk_shear(prs0/100*units.hPa, u0*units('m/s'), v0*units('m/s'), height=zh*1000*units.m, bottom=10*units.m, depth=6000*units.m)
ushr06 = shr6km[0].magnitude
vshr06 = shr6km[1].magnitude
shr06 = np.sqrt(ushr06**2 + vshr06**2)
angSHR = 180 + np.arctan2(ushr06, vshr06)*180/np.pi
shrE = mc.bulk_shear(prs0/100*units.hPa, u0*units('m/s'), v0*units('m/s'), height=zh*1000*units.m, bottom=ebot*units.m, depth=(zel-ebot)*units.m)
uEBS = shrE[0].magnitude
vEBS = shrE[1].magnitude
ebs = np.sqrt(uEBS**2 + vEBS**2)
angEBS = 180 + np.arctan2(uEBS, vEBS)*180/np.pi

srh5 = mc.storm_relative_helicity(zh*1000*units.m, u0*units('m/s'), v0*units('m/s'), depth=500*units.m, storm_u=uBR*units('m/s'), storm_v=vBR*units('m/s'))
srh500 = srh5[2].magnitude
srh1 = mc.storm_relative_helicity(zh*1000*units.m, u0*units('m/s'), v0*units('m/s'), depth=1000*units.m, storm_u=uBR*units('m/s'), storm_v=vBR*units('m/s'))
srh1km = srh1[2].magnitude
srh3 = mc.storm_relative_helicity(zh*1000*units.m, u0*units('m/s'), v0*units('m/s'), depth=3000*units.m, storm_u=uBR*units('m/s'), storm_v=vBR*units('m/s'))
srh3km = srh3[2].magnitude
srhe = mc.storm_relative_helicity(zh*1000*units.m, u0*units('m/s'), v0*units('m/s'), depth=(etop-ebot)*units.m, bottom=ebot*units.m, storm_u=uBR*units('m/s'), storm_v=vBR*units('m/s'))
esrh = srhe[2].magnitude

supcom = mc.supercell_composite(mucape*units('J/kg'), esrh*units('m^2/s^2'), ebs*units('m/s'))
scp = supcom[0].magnitude
sigtor = mc.significant_tornado(sbcape*units('J/kg'), zlcl*units.m, srh500*units('m^2/s^2'), shr06*units('m/s'))
stp = sigtor[0].magnitude



print("---BASE STATE PROFILE---")
print(f"Profile top:          {1000*max(zh)+140:.0f} m")
print(f"MUCAPE,MUCIN:         {mucape:.0f} J/kg, {mucin:.0f} J/kg")
print(f"SBCAPE,SBCIN:         {sbcape:.0f} J/kg, {sbcin:.0f} J/kg")
print(f"MLCAPE,MLCIN:         {mlcape:.0f} J/kg, {mlcin:.0f} J/kg")
print(f"DCAPE:                {dcape:.0f} J/kg")
print(f"LCL height:           {zlcl:.0f} m ({plcl:.0f} hPa)")
print(f"Eff. inflow layer:    {ebot:.0f} m-{etop:.0f} m")
print(f"Bunkers RM:           {smBR:.1f} m/s at {angBR:.0f} deg (Vector: {uBR:.1f} m/s, {vBR:.1f} m/s)")
print(f"0-6 km mean wind:     {VH06:.1f} m/s at {ang06:.0f} deg (Vector: {u06:.1f} m/s, {v06:.1f} m/s)")
print(f"-----")
print(f"Bulk shear (0-6 km):  {shr06:.1f} m/s at {angSHR:.0f} deg (Vector: {ushr06:.1f} m/s, {vshr06:.1f} m/s)")
print(f"Bulk shear (Eff.):    {ebs:.1f} m/s at {angEBS:.0f} deg (Vector: {uEBS:.1f} m/s, {vEBS:.1f} m/s)")
print(f"SRH (0-500 m):        {srh500:.0f} m2/s2")
print(f"SRH (0-1 km):         {srh1km:.0f} m2/s2")
print(f"SRH (0-3 km):         {srh3km:.0f} m2/s2")
print(f"SRH (Eff.):           {esrh:.0f} m2/s2")
print(f"Supercell composite:  {scp:.1f}")
print(f"Significant tornado:  {stp:.1f}")


fig = plt.figure(figsize=(8,8))

skew = SkewT(fig=fig)
skew.plot(prs0/100., (T0-273.15), '-r', linewidth=2)
skew.plot(prs0/100., (Td0-273.15), '-g', linewidth=2)
skew.plot(prs0/100., np.array(T_parcel.magnitude[:])-273.15, '-k', linewidth=2)
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 30)
plt.title('Base state')
ax_hod = inset_axes(skew.ax, '42%', '42%', loc=1)
H = Hodograph(ax_hod, component_range=60.)
H.add_grid(increment=20)
H.plot(u0, v0, color='k', linewidth=1.5)

plt.show()
















