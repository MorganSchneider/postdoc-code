# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 09:57:54 2025

@author: mschne28
"""

from CM1utils import *

#%% Overview plotting - dbz and thpert

fp = 'C:/Users/mschne28/Documents/cm1out/wk_p3_250m/'
# fp = 'C:/Users/mschne28/Documents/cm1r21.1/run/'
fn = np.linspace(21,28,8)


titlestr = "WK profile, dx=250m"


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
                 ds.variables['qi'][:].data[0,0,:,:]))
    thr0 = ds.variables['th0'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,:,:])
    thrpert = thr - thr0
    del thr,thr0
    
    u_gr = uinterp + ds.variables['umove'][:].data[0]
    v_gr = vinterp + ds.variables['vmove'][:].data[0]
    
    ds.close()
    
    
    xl = [-150,150]
    yl = [-150,150]
    
    xl = [-100,100]
    yl = [-100,100]
    
    
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
    ax[i,j].set_title(f"t = {time:.0f} s")
    fig.suptitle(f"Sfc dbz + max 0-2 km w ({titlestr})")
    if (n==len(fn)-1) & (figsave):
        fig.savefig(fp+f"figs/dbz.png", dpi=300)
    
    
    qix = 60
    
    plot_contourf(xh, yh, thpert, 'thpert', ax1[i,j], levels=np.linspace(-10,10,21),
                  datalims=[-10,10], xlims=xl, ylims=yl, cmap='balance', cbar=cb_flag)
    ax1[i,j].contour(xh, yh, np.max(zvort, axis=0), levels=[0.02], colors='k', linestyles='-', linewidths=1)
    ax1[i,j].quiver(xh[::qix], yh[::qix], u_gr[0,::qix,::qix], v_gr[0,::qix,::qix], color='k', scale=150, width=0.008, pivot='middle')
    ax1[i,j].set_title(f"t = {time:.0f} s")
    fig1.suptitle(f"Sfc thpert + sfc wind + max 0-1 km zeta=0.02 s$^{{-1}}$ ({titlestr})")
    if (n==len(fn)-1) & (figsave):
        fig1.savefig(fp+f"figs/thpert.png", dpi=300)


#%% Plot swaths

fp = 'C:/Users/mschne28/Documents/cm1out/wk_p3_250m/'

titlestr = 'WK profile, dx=250m'

ds = nc.Dataset(fp+"cm1out_000029.nc")
time = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
dbz = ds.variables['dbz'][:].data[0,0,:,:]
wspd = np.sqrt(ds.variables['uinterp'][:].data[0,0,:,:]**2 + ds.variables['vinterp'][:].data[0,0,:,:]**2)

sws = ds.variables['sws'][:].data[0,:,:] #max sfc wind speed
sws2 = ds.variables['sws2'][:].data[0,:,:] #translated max sfc wind speed
svs = ds.variables['svs'][:].data[0,:,:] #max sfc zeta
sps = ds.variables['sps'][:].data[0,:,:] #min sfc p'
sus = ds.variables['sus'][:].data[0,:,:] #max 5km w
shs = ds.variables['shs'][:].data[0,:,:] #max integrated UH

usfc = ds.variables['uinterp'][:].data[0,0,:,:]
vsfc = ds.variables['vinterp'][:].data[0,0,:,:]
winterp = ds.variables['winterp'][:].data[0,:,:,:]
prspert = ds.variables['prs'][:].data[0,:,:,:] - ds.variables['prs0'][:].data[0,:,:,:]
zvort = ds.variables['zvort'][:].data[0,:,:,:]
xvort = ds.variables['xvort'][:].data[0,:,:,:]
yvort = ds.variables['yvort'][:].data[0,:,:,:]
hvort = np.sqrt(xvort**2 + yvort**2)
vort = np.sqrt(xvort**2 + yvort**2 + zvort**2)
ds.close()


xl = [-150,150]
yl = [-150,150]

# xl = [-100,100]
# yl = [-100,100]

figsave = False

# print(np.max(wspd))



# fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
# plot_contourf(xh, yh, np.ma.masked_array(dbz, dbz<1), 'dbz', ax, levels=np.linspace(0,70,15),
#               datalims=[0,70], xlims=xl, ylims=yl, cmap='HomeyerRainbow', cbfs=14)
# ax.contour(xh, yh, wspd, levels=[20,26,33], colors='k', linewidths=[0.75,1,2])
# ax.set_title(f"Dbz, sfc wind speed at {time:.0f} s (NO HEAT SINK)", fontsize=16)
# # plt.show()
                      

fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, sws, 'wspd', ax, levels=np.linspace(0,40,21), 
              datalims=[0,40], xlims=xl, ylims=yl, cmap='Reds', cbfs=14)
# ax.contour(xh, yh, svs, levels=[0.01,0.02], colors=[''])
ax.contour(xh, yh, sws, levels=[26,33], colors='k', linewidths=[1,2])
ax.set_title(f"7-h sfc wind swath ({titlestr})", fontsize=16)
if figsave:
    plt.savefig(fp+"figs/wspd_swath.png", dpi=300)


# fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
# plot_contourf(xh, yh, np.ma.masked_array(svs, svs<0.001), 'zvort', ax,
#               levels=np.append(np.linspace(0,0.02,21), [0.05]), datalims=[0,0.02],
#               xlims=xl, ylims=yl, cmap='Reds', cbfs=14)
# ax.set_title(f"7-h sfc zeta swath (HEAT SINK 3-4 h)", fontsize=16)
# if figsave:
#     plt.savefig(fp+"figs/zeta_swath.png", dpi=300)


# fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
# plot_contourf(xh, yh, dbz, 'dbz', ax, levels=np.linspace(0,70,8), datalims=[0,70], 
#               xlims=xl, ylims=yl, cmap='Greys')
# ax.contour(xh, yh, sws, levels=[26,33], colors=['lightcoral','r'], linewidths=1.5)
# ax.contour(xh, yh, svs, levels=[0.01,0.02], colors=['dodgerblue','b'], linewidths=1.5)
# ax.set_title(f"dbz, max sfc wspd, max sfc zeta -- dx=1000 m", fontsize=16)


plt.show()






#%% plot cm1 stats time series

fp = 'C:/Users/mschne28/Documents/cm1out/wk_p3_250m/'

titlestr = 'WK profile, dx=250m'

ds = nc.Dataset(fp+'cm1out_000001.nc')
prs0 = ds.variables['prs0'][:].data[0,0,0,0]
ds.close()

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
sprsmax = ds.variables['sprsmax'][:].data - prs0 #max sfc p'
sprsmin = ds.variables['sprsmin'][:].data - prs0 #min sfc p'
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

l1,= ax[0].plot([0,25200], [20,20], 'tab:orange', linewidth=1.25)
l2,= ax[0].plot([0,25200], [26,26], 'red', linewidth=1.25)
l3,= ax[0].plot([0,25200], [33,33], 'darkviolet', linewidth=1.25)
ax[0].plot(time, swspmax, 'k', linewidth=2)
ax[0].set_xlim([0,25200])
ax[0].set_ylim([0,50])
# ax[0].set_xlabel('Time (s)', fontsize=14)
ax[0].set_ylabel('wspd (m/s)', fontsize=14)
ax[0].tick_params(axis='both', labelsize=12)
ax[0].set_title(f"Max sfc wind speed ({titlestr})", fontsize=16)
ax[0].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[0].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[0].yaxis.set_major_locator(MultipleLocator(10))
ax[0].yaxis.set_minor_locator(MultipleLocator(5))
ax[0].legend(handles=[l1,l2,l3], labels=['sub-severe', 'severe', 'sig severe'],
             loc='upper left', fontsize=14)

# l1,= ax[1].plot(time, sprsmax/100, 'r', linewidth=2)
# l2,= ax[1].plot(time, sprsmin/100, 'b', linewidth=2)
# ax[1].set_xlim([0,25200])
# ax[1].set_ylim([-6,8])
# ax[1].set_xlabel('Time (s)', fontsize=14)
# ax[1].set_ylabel("p' (mb)", fontsize=14)
# ax[1].tick_params(axis='both', labelsize=12)
# ax[1].set_title(f"max/min sfc prspert ({fp[35:-1]})", fontsize=16)
# ax[1].grid(visible=True, which='major', color='darkgray', linestyle='-')
# ax[1].grid(visible=True, which='minor', color='lightgray', linestyle='-')
# ax[1].xaxis.set_major_locator(MultipleLocator(3600))
# ax[1].xaxis.set_minor_locator(MultipleLocator(900))
# ax[1].yaxis.set_major_locator(MultipleLocator(2))
# ax[1].yaxis.set_minor_locator(MultipleLocator(1))
# ax[1].axhline(0, color='k', linewidth=1)

l1,= ax[1].plot(time, movmean(vortsfc,5), 'lightcoral', linewidth=2)
l2,= ax[1].plot(time, movmean(vort1km,5), 'r', linewidth=2)
l3,= ax[1].plot(time, movmean(vort3km,5), 'maroon', linewidth=2)
ax[1].set_xlim([0,25200])
ax[1].set_ylim([0,0.16])
ax[1].set_xlabel('Time (s)', fontsize=14)
ax[1].set_ylabel("zvort (1/s)", fontsize=14)
ax[1].tick_params(axis='both', labelsize=12)
ax[1].set_title(f"Max zeta ({titlestr})", fontsize=16)
ax[1].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[1].grid(visible=True, which='minor', color='lightgray', linestyle='-')
ax[1].xaxis.set_major_locator(MultipleLocator(3600))
ax[1].xaxis.set_minor_locator(MultipleLocator(900))
ax[1].yaxis.set_major_locator(MultipleLocator(0.02))
ax[1].yaxis.set_minor_locator(MultipleLocator(0.01))
ax[1].legend(handles=[l1,l2,l3], labels=['sfc','1km','3km'],
             loc='upper left', fontsize=14)

if figsave:
    plt.savefig(fp+'figs/wspd_zeta_ts.png', dpi=300)

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



#%% 

fp = 'C:/Users/mschne28/Documents/cm1out/wk_p3_250m/'

titlestr = 'WK profile, dx=250m'























