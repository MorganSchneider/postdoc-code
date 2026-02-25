# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 09:57:54 2025

@author: mschne28
"""

from CM1utils import *

#%% Overview plotting - dbz and thpert

fp = 'C:/Users/mschne28/Documents/cm1out/cwe/freeslip_wk_250m/'

fn = np.linspace(1,17,5)


if 'semislip' in fp:
    bbc = 'Semi-slip'
    sim = 'SEMISLIP'
elif 'freeslip' in fp:
    bbc = 'Free-slip'
    sim = 'FREESLIP'
elif 'noslip' in fp:
    bbc = 'No-slip'
    sim = 'NOSLIP'

titlestr = f"{bbc}, WK profile, dx=250m"
# titlestr = "New P3 -- Fir, modded"


figsave = False

fig,ax = plt.subplots(2, 4, figsize=(9.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
fig1,ax1 = plt.subplots(2, 4, figsize=(9.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')

for f in fn:
    ds = nc.Dataset(fp+f"cm1out_{f:06.0f}.nc")
    time = ds.variables['time'][:].data[0]
    xh = ds.variables['xh'][:].data
    yh = ds.variables['yh'][:].data
    zh = ds.variables['zh'][:].data
    iz1 = np.where(zh>1)[0][0]
    iz2 = np.where(zh>2)[0][0]
    iz3 = np.where(zh>3)[0][0]
    
    
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
                 ds.variables['qi1'][:].data[0,0,:,:] +
                 ds.variables['qi2'][:].data[0,0,:,:] + 
                 ds.variables['qi3'][:].data[0,0,:,:]))
                 # ds.variables['qi4'][:].data[0,0,:,:]))
    thr0 = ds.variables['th0'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,:,:])
    thrpert = thr - thr0
    del thr,thr0
    
    u_gr = uinterp + ds.variables['umove'][:].data[0]
    v_gr = vinterp + ds.variables['vmove'][:].data[0]
    
    ds.close()
    
    
    xl = [-150,150]
    yl = [-150,150]
    
    # xl = [-100,100]
    # yl = [-100,100]
    
    
    n = (f-fn[0])/(fn[1]-fn[0])
    i = int(np.floor(n/4))
    j = int(np.mod(n,4))
    
    if j == 3:
        cb_flag = True
    else:
        cb_flag = False
    
    
    plot_contourf(xh, yh, np.ma.masked_array(dbz, dbz<0.1), 'dbz', ax[i,j], levels=np.linspace(0,70,15),
                  datalims=[0,70], xlims=xl, ylims=yl, cmap='HomeyerRainbow', cbar=cb_flag)
    ax[i,j].contour(xh, yh, np.max(winterp, axis=0), levels=[5,10], colors='k', linestyles='-', linewidths=[0.5,1])
    # ax[i,j].contour(xh, yh, np.min(winterp, axis=0), levels=[-5], colors='b', linestyles='-', linewidths=1)
    if n == 0:
        l1, = ax[0,0].plot([190,200], [190,200], '-k', linewidth=0.5)
        l2, = ax[0,0].plot([190,200], [190,200], '-k', linewidth=1)
        ax[0,0].legend(handles=[l1,l2], labels=['w=5 m/s','w=10 m/s'], loc='upper right', fontsize=10)
    ax[i,j].set_title(f"t = {time:.0f} s")
    # fig.suptitle(f"Sfc dbz + max 0-2 km w ({titlestr})")
    fig.suptitle(f"{sim}")
    if (n==len(fn)-1) & (figsave):
        fig.savefig(fp+f"figs/dbz.png", dpi=300)
    
    
    qix = 60
    
    plot_contourf(xh, yh, thrpert, 'thpert', ax1[i,j], levels=np.linspace(-12,12,25),
                  datalims=[-12,12], xlims=xl, ylims=yl, cmap='balance', cbar=cb_flag)
    ax1[i,j].contour(xh, yh, np.max(zvort, axis=0), levels=[0.025], colors='r', linestyles='-', linewidths=1)
    ax1[i,j].quiver(xh[::qix], yh[::qix], u_gr[0,::qix,::qix], v_gr[0,::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
    if n == 0:
        l3, = ax1[0,0].plot([190,200], [190,200], '-r', linewidth=1)
        ax1[0,0].legend(handles=[l3], labels=["\u03B6=0.025 s$^{-1}$"], loc='upper right', fontsize=10)
    ax1[i,j].set_title(f"t = {time:.0f} s")
    # fig1.suptitle(f"Sfc thrpert + sfc wind + max 0-1 km zeta=0.025 s$^{{-1}}$ ({titlestr})")
    fig1.suptitle(f"{sim}")
    if (n==len(fn)-1) & (figsave):
        fig1.savefig(fp+f"figs/thrpert.png", dpi=300)


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

ds = nc.Dataset(fp+"cm1out_000029.nc")
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
ds.close()


xl = [-150,150]
yl = [-150,150]


figsave = False

# print(np.max(wspd))



                      

fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
c = plot_contourf(xh, yh, sws, 'wspd', ax,
              levels=np.linspace(0,40,21), datalims=[0,40], 
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,40,11), extend='max')
# ax.contour(xh, yh, svs, levels=[0.01,0.02], colors=[''])
ax.contour(xh, yh, sws, levels=[26,33], colors='k', linewidths=[0.5,1.25])
# ax.set_title(f"7-h sfc wind swath ({titlestr})", fontsize=16)
ax.set_title(f"{sim} - Surface wind swath", fontsize=16)
l1, = ax.plot([190,200], [190,200], '-k', linewidth=0.5)
l2, = ax.plot([190,200], [190,200], '-k', linewidth=1)
ax.legend(handles=[l1,l2], labels=['Severe (26 m/s)','Sig. Severe (33 m/s)'], loc='upper right', fontsize=12)
if figsave:
    plt.savefig(fp+"figs/wspd_swath.png", dpi=300)


fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, np.ma.masked_array(svs, svs<0.001), 'zvort', ax,
              levels=np.linspace(0,0.05,21), datalims=[0,0.05],
              xlims=xl, ylims=yl, cmap='Reds', cbfs=14, cbticks=np.linspace(0,0.05,11), extend='max')
# ax.set_title(f"7-h sfc zeta swath ({titlestr})", fontsize=16)
ax.set_title(f"{sim} - Surface vorticity swath", fontsize=16)
# ax.contour(xh, yh, sws, levels=[26,33], colors='k', linewidths=[0.5,1])
# l1, = ax.plot([0,0], [190,200], '-k', linewidth=0.5)
# l2, = ax.plot([0,0], [190,200], '-k', linewidth=1)
# ax.legend(handles=[l1,l2], labels=['Severe (26 m/s)','Sig. Severe (33 m/s)'], loc='upper right', fontsize=12)
ax.contour(xh, yh, sus, levels=[20,40], colors='k', linewidths=[0.75,1.5])
l1, = ax.plot([190,200], [190,200], '-k', linewidth=0.75)
l2, = ax.plot([190,200], [190,200], '-k', linewidth=1.5)
ax.legend(handles=[l1,l2], labels=['w=20 m/s','w=40 m/s'], loc='upper right', fontsize=12)
if figsave:
    plt.savefig(fp+"figs/zeta_swath.png", dpi=300)


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






#%% plot cm1 stats time series

fp = 'C:/Users/mschne28/Documents/cm1out/cwe/freeslip_wk_250m/'


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
vort3km_fs = ds.variables['vort3km'][:].data #max 3km vort
ds.close()


fp = 'C:/Users/mschne28/Documents/cm1out/cwe/noslip_wk_250m/'
ds = nc.Dataset(fp+f"cm1out_stats.nc")
wmax500_ns = ds.variables['wmax500'][:].data[:421] #max w at 500 m
wmax1000_ns = ds.variables['wmax1000'][:].data[:421] #max w at 1000 m
wmax2500_ns = ds.variables['wmax2500'][:].data[:421] #max w at 2500 m
wmax5000_ns = ds.variables['wmax5000'][:].data[:421] #max w at 5000 m
swspmax_ns = ds.variables['swspmax'][:].data[:421] #max sfc wspd
vortsfc_ns = ds.variables['vortsfc'][:].data[:421] #max sfc vort
vort1km_ns = ds.variables['vort1km'][:].data[:421] #max 1km vort
vort3km_ns = ds.variables['vort3km'][:].data[:421] #max 3km vort
ds.close()


figsave = False



fig,ax = plt.subplots(3, 1, figsize=(10,9), sharex=True, layout='constrained')

l1,= ax[0].plot(time, movmean(vortsfc_fs,5), 'k', linewidth=2)
l2,= ax[0].plot(time, movmean(vortsfc_ns,5), 'dodgerblue', linewidth=2)
l3,= ax[0].plot(time, movmean(vortsfc,5), 'crimson', linewidth=2)
# l1,= ax[0].plot(time, vortsfc_fs, 'k', linewidth=2)
# l2,= ax[0].plot(time,vortsfc_ns, 'dodgerblue', linewidth=2)
# l3,= ax[0].plot(time, vortsfc, 'crimson', linewidth=2)
ax[0].set_xlim([0,25200])
ax[0].set_ylim([0,0.3])
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
ax[1].set_ylim([0,0.25])
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
    plt.savefig('C:/Users/mschne28/Documents/cm1out/cwe/zeta_all_timeseries.png', dpi=300)
# plt.show()

#%


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
ax[2].set_ylim([20,50])
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
    plt.savefig('C:/Users/mschne28/Documents/cm1out/cwe/w_all_timeseries.png', dpi=300)

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





