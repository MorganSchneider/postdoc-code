# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 09:01:52 2025

@author: mschne28
"""

from CM1utils import *

#%% plot dbz, thpert, w, zeta, etc

# semi-slip test1 and test2 are the most interesting/bow echo-y
# Cd = 0.0014 and 0.005
# fp = 'C:/Users/mschne28/Documents/cm1out/semislip_wk_250m/'
fp = 'C:/Users/mschne28/Documents/cm1r21.1/run/'
fn = np.linspace(2,9,8)

p3_version = 'old'

# titlestr = "Semi-slip, dx=250m"
titlestr = f"Test - {p3_version} P3"


figsave = False

fig,ax = plt.subplots(2, 4, figsize=(9.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
fig1,ax1 = plt.subplots(2, 4, figsize=(9.5,5), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')

for f in fn:
    ds = nc.Dataset(fp+f"cm1out_{f:06.0f}_{p3_version}P3.nc")
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
    
    # thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] - 
    #             (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] + 
    #              ds.variables['qi'][:].data[0,0,:,:] + ds.variables['qs'][:].data[0,0,:,:] + 
    #              ds.variables['qg'][:].data[0,0,:,:] + ds.variables['qhl'][:].data[0,0,:,:]))
    # thr0 = ds.variables['th0'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,:,:])
    # thrpert = thr - thr0
    # del thr,thr0
    
    u_gr = uinterp + ds.variables['umove'][:].data[0]
    v_gr = vinterp + ds.variables['vmove'][:].data[0]
    
    ds.close()
    
    
    xl = [-150,150]
    yl = [-150,150]
    
    xl = [-60,60]
    yl = [-60,60]
    
    
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
    
    
    qix = 3
    
    plot_contourf(xh, yh, thpert, 'thpert', ax1[i,j], levels=np.linspace(-8,8,17), #was levels=np.linspace(-15,15,31), datalims=[-15,15]
                  datalims=[-8,8], xlims=xl, ylims=yl, cmap='balance', cbar=cb_flag)
    ax1[i,j].contour(xh, yh, np.max(zvort, axis=0), levels=[0.02], colors='k', linestyles='-', linewidths=1)
    ax1[i,j].quiver(xh[::qix], yh[::qix], u_gr[0,::qix,::qix], v_gr[0,::qix,::qix], color='k', scale=150, width=0.005, pivot='middle')
    ax1[i,j].set_title(f"t = {time:.0f} s")
    fig1.suptitle(f"Sfc thpert + sfc winds + max 0-1 km zeta=0.02 s$^{{-1}}$ ({titlestr})")
    if (n==len(fn)-1) & (figsave):
        fig1.savefig(fp+f"figs/thpert.png", dpi=300)




#%% Plot tendency bubbles over dbz


fp = 'C:/Users/mschne28/Documents/cm1out/semislip_500m/'

ds = nc.Dataset(fp+'cm1out_000013.nc')
time1 = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['zh'][:].data
iz = np.argmin(abs(zh-1))
dbz1 = ds.variables['dbz'][:].data[0,0,:,:]
thpert1 = ds.variables['th'][:].data[0,0,:,:] - ds.variables['th0'][:].data[0,0,:,:]
ugr1 = ds.variables['uinterp'][:].data[0,0,:,:] + ds.variables['umove'][:].data[0]
vgr1 = ds.variables['vinterp'][:].data[0,0,:,:] + ds.variables['vmove'][:].data[0]
winterp1 = ds.variables['winterp'][:].data[0,0:iz,:,:]
zvort1 = ds.variables['zvort'][:].data[0,0:iz,:,:]
ds.close()

ds = nc.Dataset(fp+'cm1out_000017.nc')
time2 = ds.variables['time'][:].data[0]
dbz2 = ds.variables['dbz'][:].data[0,0,:,:]
thpert2 = ds.variables['th'][:].data[0,0,:,:] - ds.variables['th0'][:].data[0,0,:,:]
ugr2 = ds.variables['uinterp'][:].data[0,0,:,:] + ds.variables['umove'][:].data[0]
vgr2 = ds.variables['vinterp'][:].data[0,0,:,:] + ds.variables['vmove'][:].data[0]
winterp2 = ds.variables['winterp'][:].data[0,0:iz,:,:]
zvort2 = ds.variables['zvort'][:].data[0,0:iz,:,:]
ds.close()






bigR = np.zeros(shape=(len(yh),len(xh)), dtype=float)
coldxc = -40
coldyc = 35
coldrad = 30
coldrady = 25
s0 = -0.010

xrad = np.abs(xh-coldxc)
yrad = np.abs(yh-coldyc)

xi = (xrad < coldrad)

for j in range(len(yh)):
    if yrad[j] < coldrady:
        bigR[j,xi] = 1 - (xrad[xi]**2 / coldrad**2)
    elif (yrad[j] >= coldrady) & (yrad[j] < coldrad+coldrady):
        bigR[j,xi] = 1 - ( (xrad[xi]**2 + (yrad[j]-coldrady)**2) / coldrad**2)


s = s0*bigR
slevs = [0.001, 0.25, 0.5, 0.8, 0.99]

xl = [-150,150]
yl = [-150,150]





fig,ax = plt.subplots(1, 2, figsize=(10,4), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')

plot_contourf(xh, yh, np.ma.masked_array(dbz1, dbz1<1), 'dbz', ax[0],
              levels=np.linspace(0,70,15), datalims=[0,70], xlims=xl, ylims=yl,
              cmap='HomeyerRainbow', cbar=False)
ax[0].contour(xh, yh, np.max(winterp1, axis=0), levels=[5,10], colors='k', linestyles='-', linewidths=1)
ax[0].contour(xh, yh, bigR, levels=slevs, colors='k', linewidths=0.75)
ax[0].set_title(f"t = {time1:.0f} s", fontsize=14)
ax[0].grid(visible=True, which='both')
ax[0].xaxis.set_major_locator(MultipleLocator(50))
ax[0].xaxis.set_minor_locator(MultipleLocator(25))
ax[0].yaxis.set_major_locator(MultipleLocator(50))
ax[0].yaxis.set_minor_locator(MultipleLocator(25))

plot_contourf(xh, yh, np.ma.masked_array(dbz2, dbz2<1), 'dbz', ax[1],
              levels=np.linspace(0,70,15), datalims=[0,70], xlims=xl, ylims=yl,
              cmap='HomeyerRainbow', cbar=True)
ax[1].contour(xh, yh, np.max(winterp2, axis=0), levels=[5,10], colors='k', linestyles='-', linewidths=1)
ax[1].contour(xh, yh, bigR, levels=slevs, colors='k', linewidths=0.75)
ax[1].set_title(f"t = {time2:.0f} s", fontsize=14)
ax[1].grid(visible=True, which='both')
ax[1].xaxis.set_major_locator(MultipleLocator(50))
ax[1].xaxis.set_minor_locator(MultipleLocator(25))
ax[1].yaxis.set_major_locator(MultipleLocator(50))
ax[1].yaxis.set_minor_locator(MultipleLocator(25))
plt.show()


qix = 30


fig,ax = plt.subplots(1, 2, figsize=(10,4), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')

plot_contourf(xh, yh, thpert1, 'thpert', ax[0],
              levels=np.linspace(-10,10,21), datalims=[-10,10], xlims=xl, ylims=yl,
              cmap='balance', cbar=False)
ax[0].contour(xh, yh, np.max(zvort1, axis=0), levels=[0.01], colors='k', linestyles='-', linewidths=1)
ax[0].contour(xh, yh, bigR, levels=slevs, colors='k', linewidths=0.75)
ax[0].quiver(xh[::qix], yh[::qix], ugr1[::qix,::qix], vgr1[::qix,::qix], color='k', scale=200, width=0.005, pivot='middle')
ax[0].set_title(f"t = {time1:.0f} s", fontsize=14)
ax[0].grid(visible=True, which='both')
ax[0].xaxis.set_major_locator(MultipleLocator(50))
ax[0].xaxis.set_minor_locator(MultipleLocator(25))
ax[0].yaxis.set_major_locator(MultipleLocator(50))
ax[0].yaxis.set_minor_locator(MultipleLocator(25))

plot_contourf(xh, yh, thpert2, 'thpert', ax[1],
              levels=np.linspace(-10,10,21), datalims=[-10,10], xlims=xl, ylims=yl,
              cmap='balance', cbar=True)
ax[1].contour(xh, yh, np.max(zvort2, axis=0), levels=[0.01], colors='k', linestyles='-', linewidths=1)
ax[1].contour(xh, yh, bigR, levels=slevs, colors='k', linewidths=0.75)
ax[1].quiver(xh[::qix], yh[::qix], ugr2[::qix,::qix], vgr2[::qix,::qix], color='k', scale=200, width=0.005, pivot='middle')
ax[1].set_title(f"t = {time2:.0f} s", fontsize=14)
ax[1].grid(visible=True, which='both')
ax[1].xaxis.set_major_locator(MultipleLocator(50))
ax[1].xaxis.set_minor_locator(MultipleLocator(25))
ax[1].yaxis.set_major_locator(MultipleLocator(50))
ax[1].yaxis.set_minor_locator(MultipleLocator(25))
plt.show()


#%% cross sections

fp = 'C:/Users/mschne28/Documents/cm1out/semislip_500m/'

ds = nc.Dataset(fp+'cm1out_000023.nc')
time = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['zh'][:].data
iz = np.where(zh >= 10)[0][1]
yc = -85
iy = np.argmin(abs(yh-yc))
xc = -18
ix = np.argmin(abs(xh-xc))

dbz = ds.variables['dbz'][:].data[0,0:iz,:,:]
thpert = ds.variables['th'][:].data[0,0:iz,:,:] - ds.variables['th0'][:].data[0,0:iz,:,:]
uinterp = ds.variables['uinterp'][:].data[0,0:iz,:,:]
vinterp = ds.variables['vinterp'][:].data[0,0:iz,:,:]
winterp = ds.variables['winterp'][:].data[0,0:iz,:,:]
zvort = ds.variables['zvort'][:].data[0,0:iz,:,:]
xvort = ds.variables['xvort'][:].data[0,0:iz,:,:]
yvort = ds.variables['yvort'][:].data[0,0:iz,:,:]

# thr = ds.variables['th'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv'][:].data[0,0,:,:] - 
#             (ds.variables['qc'][:].data[0,0,:,:] + ds.variables['qr'][:].data[0,0,:,:] + 
#              ds.variables['qi'][:].data[0,0,:,:] + ds.variables['qs'][:].data[0,0,:,:] + 
#              ds.variables['qg'][:].data[0,0,:,:] + ds.variables['qhl'][:].data[0,0,:,:]))
# thr0 = ds.variables['th0'][:].data[0,0,:,:] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,:,:])
# thrpert = thr - thr0
# del thr,thr0

u_gr = uinterp + ds.variables['umove'][:].data[0]
v_gr = vinterp + ds.variables['vmove'][:].data[0]

ds.close()


xl = [-50,0]
yl = [-100,-50]
 

fig,ax = plt.subplots(1, 1, figsize=(8,6), layout='constrained')
plot_contourf(xh, yh, thpert[0,:,:], 'thpert', ax, levels=np.linspace(-12,12,25),
              datalims=[-12,12], xlims=xl, ylims=yl, cmap='balance')
ax.axhline(yc, color='k', linewidth=1.5)
ax.axvline(xc, color='k', linewidth=1.5)
ax.quiver(xh[::4], yh[::4], u_gr[0,::4,::4], v_gr[0,::4,::4], color='k', scale=300, width=0.0035, pivot='middle')
ax.set_title(f"t = {time:.0f} s")
plt.show()

# if figsave:
#     plt.savefig(fp+f"figs/thpert_{time:.0f}s.png", dpi=300)



fig,ax = plt.subplots(1, 1, figsize=(9,5), layout='constrained')
# plot_contourf(xh, zh[:iz], np.ma.masked_array(dbz[:,iy,:], dbz[:,iy,:]<1), 'dbz', ax, levels=np.linspace(0,70,15),
#               datalims=[0,70], xlims=[-100,100], ylims=[0,10], cmap='HomeyerRainbow')
# plot_contourf(xh, zh[:iz], thpert[:,iy,:], 'thpert', ax, levels=np.linspace(-12,12,15),
#               datalims=[-12,12], xlims=[-50,50], ylims=[0,5], cmap='balance')
plot_contourf(xh, zh[:iz], winterp[:,iy,:], 'w', ax, levels=np.linspace(-30,30,13),
              datalims=[-30,30], xlims=[-30,0], ylims=[0,5], cmap='balance')
ax.quiver(xh[::3], zh[:iz][::4], u_gr[::4,iy,::3], winterp[::4,iy,::3], color='k', scale=400, width=0.0035, pivot='middle')
ax.contour(xh, zh[:iz], thpert[:,iy,:], levels=[-8,-5], colors='k', linestyles='-', linewidths=1.5)
plt.show()


fig,ax = plt.subplots(1, 1, figsize=(9,5), layout='constrained')
# plot_contourf(xh, zh[:iz], np.ma.masked_array(dbz[:,iy,:], dbz[:,iy,:]<1), 'dbz', ax, levels=np.linspace(0,70,15),
#               datalims=[0,70], xlims=[-100,100], ylims=[0,10], cmap='HomeyerRainbow')
# plot_contourf(yh, zh[:iz], thpert[:,:,ix], 'thpert', ax, levels=np.linspace(-12,12,15),
#               datalims=[-12,12], xlims=[-100,0], ylims=[0,5], cmap='balance')
plot_contourf(yh, zh[:iz], winterp[:,:,ix], 'w', ax, levels=np.linspace(-30,30,13),
              datalims=[-30,30], xlims=[-100,-70], ylims=[0,5], cmap='balance')
ax.quiver(yh[::3], zh[:iz][::4], v_gr[::4,::3,ix], winterp[::4,::3,ix], color='k', scale=400, width=0.0035, pivot='middle')
# ax.contour(yh, zh[:iz], xvort[:,:,ix], levels=[-0.1,-0.05,0.05,0.1], colors='k', linestyles=['--','--','-','-'], linewidths=1.5)
ax.contour(yh, zh[:iz], thpert[:,:,ix], levels=[-8,-5], colors='k', linestyles='-', linewidths=1.5)
plt.show()





# fig,ax = plt.subplots(1, 1, figsize=(9,5), layout='constrained')
# plot_contourf(xh, zh[:iz], yvort[:,iy,:], 'yvort', ax, levels=np.linspace(-0.2,0.2,41),
#               datalims=[-0.2,0.2], xlims=[-50,20], ylims=[0,2], cmap='balance')
# ax.quiver(xh[::6], zh[:iz][::3], u_gr[::3,iy,::6], winterp[::3,iy,::6], color='k', scale=500, width=0.0035, pivot='middle')
# ax.contour(xh, zh[:iz], thpert[:,iy,:], levels=[-8,-5], colors='k', linestyles='-', linewidths=1.5)
# plt.show()


fig,ax = plt.subplots(1, 1, figsize=(9,5), layout='constrained')
plot_contourf(yh, zh[:iz], zvort[:,:,ix], 'zvort', ax, levels=np.linspace(-0.05,0.05,21),
              datalims=[-0.05,0.05], xlims=[-100,-70], ylims=[0,5], cmap='balance')
ax.quiver(yh[::3], zh[:iz][::4], v_gr[::4,::3,ix], winterp[::4,::3,ix], color='k', scale=500, width=0.0035, pivot='middle')
ax.contour(yh, zh[:iz], thpert[:,:,ix], levels=[-8,-5], colors='k', linestyles='-', linewidths=1.5)
plt.show()




# fig,ax = plt.subplots(1, 1, figsize=(9,5), layout='constrained')
# plot_contourf(xh, zh[:iz], yvort[:,iy,:], 'yvort', ax, levels=np.linspace(-0.2,0.2,41),
#               datalims=[-0.2,0.2], xlims=[-15,0], ylims=[0,1], cmap='balance')
# ax.quiver(xh, zh[:iz][::2], u_gr[::2,iy,:], winterp[::2,iy,:], color='k', scale=500, width=0.0035, pivot='middle')
# ax.contour(xh, zh[:iz], thpert[:,iy,:], levels=[-8,-5], colors='k', linestyles='-', linewidths=1.5)
# plt.show()

# fig,ax = plt.subplots(1, 1, figsize=(9,5), layout='constrained')
# plot_contourf(yh, zh[:iz], xvort[:,:,ix], 'xvort', ax, levels=np.linspace(-0.2,0.2,41),
#               datalims=[-0.2,0.2], xlims=[-100,-85], ylims=[0,1], cmap='balance')
# ax.quiver(yh, zh[:iz][::2], v_gr[::2,:,ix], winterp[::2,:,ix], color='k', scale=500, width=0.0035, pivot='middle')
# ax.contour(yh, zh[:iz], thpert[:,:,ix], levels=[-8,-5], colors='k', linestyles='-', linewidths=1.5)
# plt.show()



#%% plot cm1out swaths

fp = 'C:/Users/mschne28/Documents/cm1out/warmbub_p3_250m/'

titlestr = 'P3 3-moment 250-m'

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

# xl = [-30,50]
# yl = [-150,-70]

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

fp = 'C:/Users/mschne28/Documents/cm1out/warmbub_p3_250m/'

titlestr = 'P3 3-moment 250-m'

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


#%% Testing vorticity budgets?

# something is fucked with the budget calculations and idk how to fix it


ds = nc.Dataset('C:/Users/mschne28/Documents/cm1r21.1/run/cm1out_prclvort.nc')

component = 'x'

ptime = ds.variables['time'][:].data
pid = ds.variables['xh'][:].data

vort = ds.variables[component+'vort'][:].data[:,:]
tilt = ds.variables[component+'vtilt'][:].data[:,:]
stretch = ds.variables[component+'vstretch'][:].data[:,:]
bcl = ds.variables[component+'vbcl'][:].data[:,:]
fric = ds.variables[component+'vfric'][:].data[:,:]
ten = ds.variables[component+'vten'][:].data[:,:]

ds.close()



# ten = ten/5
# fric = ten - tilt - stretch - bcl



fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.plot(ptime, np.median(tilt, axis=1), 'r')
ax.plot(ptime, np.median(stretch, axis=1), 'b')
ax.plot(ptime, np.median(bcl, axis=1), 'g')
ax.plot(ptime, np.median(fric, axis=1), 'm')
ax.plot(ptime, np.median(ten, axis=1), 'k')

plt.show()




fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.plot(ptime, tilt[:,1200], 'r')
ax.plot(ptime, stretch[:,1200], 'b')
ax.plot(ptime, bcl[:,1200], 'g')
ax.plot(ptime, fric[:,1200], 'm')
ax.plot(ptime, ten[:,1200], 'k')

plt.show()



fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.plot(ptime, np.cumsum(np.median(tilt*5, axis=1)), 'r')
ax.plot(ptime, np.cumsum(np.median(stretch*5, axis=1)), 'b')
ax.plot(ptime, np.cumsum(np.median(bcl*5, axis=1)), 'g')
ax.plot(ptime, np.cumsum(np.median(fric*5, axis=1)), 'm')
ax.plot(ptime, np.cumsum(np.median(ten*5, axis=1)), 'k')
ax.plot(ptime, np.median(vort-vort[0,:], axis=1), '--k')

plt.show()




fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.plot(ptime, np.cumsum(tilt[:,1200]*5), 'r')
ax.plot(ptime, np.cumsum(stretch[:,1200]*5), 'b')
ax.plot(ptime, np.cumsum(bcl[:,1200]*5), 'g')
ax.plot(ptime, np.cumsum(fric[:,1200]*5), 'm')
ax.plot(ptime, np.cumsum(ten[:,1200]*5), 'k')
ax.plot(ptime, vort[:,1200]-vort[0,1200], '--k')

plt.show()






#%% calculate my own swaths

sim = 'test7'

fp = f"C:/Users/mschne28/Documents/cm1out/{sim}/"

ds = nc.Dataset(fp+f"cm1out_000001.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['zh'][:].data
ds.close()

iz500 = np.argmin(abs(zh-0.5))
iz1000 = np.argmin(abs(zh-1))
iz1500 = np.argmin(abs(zh-1.5))
iz2000 = np.argmin(abs(zh-2))
iz2500 = np.argmin(abs(zh-2.5))
iz3000 = np.argmin(abs(zh-3))
iz4000 = np.argmin(abs(zh-4))
iz5000 = np.argmin(abs(zh-5))

zvort_swath = np.zeros(shape=(len(zh),len(yh),len(xh)), dtype=float)
for fn in np.arange(1,30):
    ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    zvort = ds.variables['zvort'][:].data[0,0,:,:]
    ds.close()
    
    zvort_swath = np.maximum(zvort_swath, zvort)



xl = [-150,150]
yl = [-150,150]

fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, zvort_swath[iz500,:,:], 'zvort', ax, levels=np.linspace(0,0.015,16), 
              datalims=[0,0.015], xlims=xl, ylims=yl, cmap='LangRainbow12')
ax.set_title(f"max 0.5km zeta ({sim})", fontsize=16)


plt.show()


#%% making my own long warm bubble

beta = np.zeros(shape=(len(yh),len(xh)), dtype=float)

ric = 0
rjc = 0
bhrad = 10
bubdy = 50

for i in range(len(xh)):
    for j in range(len(yh)):
        if np.abs(yh[j]-rjc) < bubdy:
            beta[j,i] = np.sqrt( ((xh[i]-ric)/bhrad)**2)
        else:
            beta[j,i] = np.sqrt(((xh[i]-ric)/bhrad)**2 + 
                           ((np.abs(yh[j]-rjc)-bubdy)/bhrad)**2)
            
fig,ax = plt.subplots(1, 1, figsize=(8,6), layout='constrained')
plot_contourf(xh, yh, 1-beta, 'dbz', ax, levels=np.linspace(0,1,11),
              datalims=[0,1], xlims=[-100,100], ylims=[-100,100], cmap='Reds')
# ax.set_title(f"t = {time:.0f} s")
plt.show()



bigR = np.zeros(shape=(len(yh),len(xh)), dtype=float)


coldxc = 0
coldyc = 0
coldzc = 0
coldrad = 10
coldrady = 40
colddepth = 3
coldamp = -0.020


for j in range(len(yh)):
    for i in range(len(xh)):
        if np.abs(xh[i]-coldxc) < coldrad:
            if np.abs(yh[j]-coldyc) < coldrady:
                radius = xh[i]-coldxc
                bigR[j,i] = 1 - (radius**2 / coldrad**2)
            elif (np.abs(yh[j]-coldyc) >= coldrady) & (np.abs(yh[j]-coldyc) < coldrady+coldrad):
                radius = np.sqrt((xh[i]-coldxc)**2 + (np.abs(yh[j]-coldyc)-coldrady)**2)
                bigR[j,i] = 1 - (radius**2 / coldrad**2)
            else:
                bigR[j,i] = 0
        else:
            bigR[j,i] = 0



fig,ax = plt.subplots(1, 1, figsize=(8,6), layout='constrained')
plot_contourf(xh, yh, bigR, 'dbz', ax, levels=np.linspace(0,1,51),
              datalims=[0,1], xlims=[-150,150], ylims=[-150,150], cmap='Reds')
# ax.set_title(f"t = {time:.0f} s")
plt.show()






#%% Weisman Klemp analytic sounding

def rslf(p,t):
    esl = 611.2*np.exp(17.67*(t-273.15)/(t-29.65))
    esl = np.minimum(esl, p*0.5)
    rslf = eps*esl/(p-esl)
    return rslf

ds = nc.Dataset('C:/Users/mschne28/Documents/cm1out/semislip_test1/cm1out_000001.nc')
zh = ds.variables['zh'][:].data*1000
zf = ds.variables['zf'][:].data*1000
ds.close()

z_trop = 12000
th_trop = 343
t_trop = 213
th_sfc = 300
prs_sfc = 100000
qv_pbl = 0.014

g = 9.81
rd = 287.04
rv = 461.5
lv = 2.501e6
cp = 1005.7
cv = cp-rd
p00 = 1e5
eps = rd/rv
reps = rv/rd

pi_sfc = (prs_sfc/p00)**(rd/cp)
qv_sfc = rslf(prs_sfc, th_sfc*pi_sfc)
thv_sfc = th_sfc*(1+qv_sfc*reps)/(1+qv_sfc)

rh0 = np.zeros(shape=(len(zh),), dtype=float)
th0 = np.zeros(shape=(len(zh),), dtype=float)
qv0 = np.zeros(shape=(len(zh),), dtype=float)
thv0 = np.zeros(shape=(len(zh),), dtype=float)
pi0 = np.zeros(shape=(len(zh),), dtype=float)
prs0 = np.zeros(shape=(len(zh),), dtype=float)
td0 = np.zeros(shape=(len(zh),), dtype=float)


for k in range(len(zh)):
    if zh[k] < z_trop:
        th0[k] = th_sfc + (th_trop-th_sfc)*(zh[k]/z_trop)**1.25
        rh0[k] = 1 - 0.75*(zh[k]/z_trop)**1.25
    else:
        th0[k] = th_trop*np.exp(g/(t_trop*cp)*(zh[k]-z_trop))
        rh0[k] = 0.25


for n in range(20):
    for k in range(len(zh)):
        thv0[k] = th0[k] * (1+reps*qv0[k])/(1+qv0[k])
    
    pi0[0] = pi_sfc - g*zh[0]/(0.5*cp*(thv_sfc+thv0[0]))
    for k in np.arange(1,len(zh)):
        pi0[k] = pi0[k-1] - g*(zh[k]-zh[k-1]) / (0.5*cp*(thv0[k]+thv0[k-1]))
    
    for k in range(len(zh)):
        prs0[k] = p00*(pi0[k]**(rd/cp))
    
    for k in range(len(zh)):
        qv0[k] = rh0[k] * rslf(prs0[k], th0[k]*pi0[k])
        if qv0[k] > qv_pbl:
            qv0[k] = qv_pbl

for k in range(len(zh)):
    rh0[k] = qv0[k] / (rslf(prs0[k], th0[k]*pi0[k]))

for k in range(len(zh)):
    es = 611.2 * np.exp(17.67*(th0[k]*pi0[k]-273.15)/(th0[k]*pi0[k]-29.65))
    e = qv0[k]*(prs0[k]-es)*reps
    td0[k] = ((1/273.15) - (rv/lv)*np.log(e/611.2))**-1

psurf = prs_sfc
tsurf = th_sfc * (psurf/p00)**(rd/cp)
qsurf = qv_pbl



udep1 = 2000
udep2 = 6000
umax1 = 7
umax2 = 31
vmax1 = umax1

u0 = np.zeros(shape=(len(zh),), dtype=float)
v0 = np.zeros(shape=(len(zh),), dtype=float)
u1 = np.zeros(shape=(len(zh),), dtype=float)
v1 = np.zeros(shape=(len(zh),), dtype=float)



for k in range(len(zh)):
    zu = zh[k]
    if zu <= udep1:
        ang = 90*(zu/udep1)*(np.pi/180)
        u0[k] = umax1 - umax1*np.cos(ang)
        v0[k] = vmax1*np.sin(ang)
        # u0[k] = umax1 - umax1*np.cos(ang) - umax1*np.sin(ang-20*np.pi/180)
    elif (zu > udep1) & (zu <= udep2):
        u0[k] = umax1 + (zu-udep1)*(umax2-umax1)/(udep2-udep1)
        v0[k] = vmax1
    else:
        u0[k] = umax2
        v0[k] = vmax1


V0 = np.sqrt(u0**2 + v0**2)

# Rotate hodograph by 30 degrees (or whatever angle)
rot_ang = -20
u1 = u0*np.cos(rot_ang*np.pi/180) + v0*np.sin(rot_ang*np.pi/180)
v1 = v0*np.cos(rot_ang*np.pi/180) - u0*np.sin(rot_ang*np.pi/180)
V1 = np.sqrt(u1**2 + v1**2)



u_06 = np.mean(u0[0:81])
v_06 = np.mean(v0[0:81])
u_shr = np.mean(u0[77:81]) - np.mean(u0[0:25])
v_shr = np.mean(v0[77:81]) - np.mean(v0[0:25])
V_shr = np.sqrt(u_shr**2 + v_shr**2)
D = 7.5

u_rm = u_06 + D*(v_shr/V_shr)
v_rm = v_06 - D*(u_shr/V_shr)

u1_06 = np.mean(u1[0:81])
v1_06 = np.mean(v1[0:81])
u1_shr = np.mean(u1[77:81]) - np.mean(u1[0:25])
v1_shr = np.mean(v1[77:81]) - np.mean(v1[0:25])
V1_shr = np.sqrt(u1_shr**2 + v1_shr**2)

u1_rm = u1_06 + D*(v1_shr/V1_shr)
v1_rm = v1_06 - D*(u1_shr/V1_shr)


# fig,ax = plt.subplots(1, 2, figsize=(6,6), sharey=True, layout='constrained')

# ax[0].plot(th0*pi0-273.15, zh/1000, 'r', linewidth=2)
# ax[0].plot(td0-273.15, zh/1000, 'green', linewidth=2)
# ax[0].grid(visible=True, which='both')
# # ax[0].set_xlim([])

# ax[1].quiver(0*zh[::2], zh[::2]/1000, u0[::2], v0[::2], color='k', scale=100, width=0.01, pivot='tail')
# ax[1].grid(visible=True, which='both')
# # ax[1].set_xlim([-2,2])
# # ax[1].set_ylim([0,6])


fig,ax = plt.subplots(3, 1, figsize=(5,5), layout='constrained')
ax[0].plot(zh, u0, 'k')
ax[0].plot(zh, u1, 'b')
ax[1].plot(zh, v0, 'k')
ax[1].plot(zh, v1, 'b')
ax[2].plot(zh, V0, 'k')
ax[2].plot(zh, V1, 'b')


phi = np.linspace(0,360,361)


fig,ax = plt.subplots(1, 1, figsize=(6,6), layout='constrained')
ax.plot(10*np.cos(phi*np.pi/180), 10*np.sin(phi*np.pi/180), 'grey', linewidth=1)
ax.plot(20*np.cos(phi*np.pi/180), 20*np.sin(phi*np.pi/180), 'grey', linewidth=1)
ax.plot(30*np.cos(phi*np.pi/180), 30*np.sin(phi*np.pi/180), 'grey', linewidth=1)
ax.plot(40*np.cos(phi*np.pi/180), 40*np.sin(phi*np.pi/180), 'grey', linewidth=1)
# ax.plot(50*np.cos(phi*np.pi/180), 50*np.sin(phi*np.pi/180), 'grey', linewidth=1)
ax.plot(u0, v0, 'k', linewidth=2)
ax.plot([u0[0],u_06], [v0[0],v_06], 'k', linewidth=3)
ax.plot([u0[0],u_shr], [v0[0],v_shr], 'b', linewidth=3)
ax.plot([u0[0],u_rm], [v0[0],v_rm], 'r', linewidth=3)
ax.plot(u1, v1, 'k', linestyle='--', linewidth=2)
ax.plot([u1[0],u1_06], [v1[0],v1_06], 'k', linestyle='--', linewidth=3)
ax.plot([u1[0],u1_shr], [v1[0],v1_shr], 'b', linestyle='--', linewidth=3)
ax.plot([u1[0],u1_rm], [v1[0],v1_rm], 'r', linestyle='--', linewidth=3)
ax.set_xlim([-30,30])
ax.set_ylim([-30,30])
ax.axhline(0, color='grey', linewidth=1)
ax.axvline(0, color='grey', linewidth=1)
plt.show()


#%% plot base state changes throughout simulations

fp = 'C:/Users/mschne28/Documents/cm1out/semislip_test1/'
fn = np.linspace(1,29,8)

fig0,ax0 = plt.subplots(2, 4, figsize=(8,7), sharex=True, sharey=True, layout='constrained')
fig1,ax1 = plt.subplots(2, 4, figsize=(8,4), sharex=True, sharey=True, layout='constrained')



for f in fn:
    ds = nc.Dataset(fp+f"cm1out_{f:06.0f}.nc")
    time = ds.variables['time'][:].data[0]
    zh = ds.variables['zh'][:].data
    i6 = np.argmin(abs(zh-6))
    
    th = ds.variables['th'][:].data[0,:,1,-2]
    qv = ds.variables['qv'][:].data[0,:,1,-2]
    prs = ds.variables['prs'][:].data[0,:,1,-2]
    t = th * (prs/100000)**(287.04/1005.7)
    e = qv*prs*461.5/287.04
    td = 206.4/(1-0.038*np.log(e))
    u = ds.variables['uinterp'][:].data[0,:i6+1,1,-2] + ds.variables['umove'][:].data[0]
    v = ds.variables['vinterp'][:].data[0,:i6+1,1,-2]
    ds.close()
    
    
    phi = np.linspace(0,360,361)
    
    n = (f-1)/4
    if n < 4:
        i = 0
        j = int(n)
    else:
        i = 1
        j = int(n-4)
    
    ax0[i,j].plot(t-273.15, zh, 'r', linewidth=2)
    ax0[i,j].plot(td-273.15, zh, 'green', linewidth=2)
    ax0[i,j].grid(visible=True, which='both')
    ax0[i,j].set_title(f"{time:.0f} s")
    
    ax1[i,j].plot(10*np.cos(phi*np.pi/180), 10*np.sin(phi*np.pi/180), 'lightgrey', linewidth=1)
    ax1[i,j].plot(20*np.cos(phi*np.pi/180), 20*np.sin(phi*np.pi/180), 'lightgrey', linewidth=1)
    ax1[i,j].plot(30*np.cos(phi*np.pi/180), 30*np.sin(phi*np.pi/180), 'lightgrey', linewidth=1)
    ax1[i,j].plot(40*np.cos(phi*np.pi/180), 40*np.sin(phi*np.pi/180), 'lightgrey', linewidth=1)
    ax1[i,j].plot(50*np.cos(phi*np.pi/180), 50*np.sin(phi*np.pi/180), 'lightgrey', linewidth=1)
    ax1[i,j].plot(u, v, 'k', linewidth=1.5)
    ax1[i,j].axhline(0, color='grey', linewidth=1)
    ax1[i,j].axvline(0, color='grey', linewidth=1)
    ax1[i,j].set_title(f"{time:.0f} s", fontsize=10)

# ax0[0,0].set_xlim([-70,30])
# ax0[0,0].set_ylim([0,15])
# ax0[0,0].text(-50, 14, f"{fp[35:-1]}", fontsize=15, color='k')
ax0[0,0].set_xlim([0,30])
ax0[0,0].set_ylim([0,3])
ax0[0,0].text(1, 2.8, f"{fp[35:-1]}", fontsize=15, color='k')
ax0[0,0].xaxis.set_major_locator(MultipleLocator(10))

# ax1[0,0].set_xlim([-40,40])
# ax1[0,0].set_ylim([-40,40])
# ax1[0,0].text(-38, 32, f"{fp[35:-1]}", fontsize=11, color='k')
ax1[0,0].set_xlim([-20,40])
ax1[0,0].set_ylim([-20,40])
ax1[0,0].text(-19, 34, f"{fp[35:-1]}", fontsize=11, color='k')
ax1[0,0].xaxis.set_major_locator(MultipleLocator(20))

if 'semislip' in fp:
    if fp[-2] == '1':
        Cd = 0.0014
    elif fp[-2] == '2':
        Cd = 0.005
    elif fp[-2] == '3':
        Cd = 0.01
    elif fp[-2] == '4':
        Cd = 0.02
    elif fp[-2] == '5':
        Cd = 0.05
    ax0[0,0].text(1, 2.53, f"Cd={Cd}", fontsize=15, color='k')
    # ax1[0,0].text(-38, 22, f"Cd={Cd}", fontsize=11, color='k')
    ax1[0,0].text(-19, 27, f"Cd={Cd}", fontsize=11, color='k')

plt.show()

#%% hodograph rotation (not useful)


ds = nc.Dataset('C:/Users/mschne28/Documents/cm1out/semislip_test7/cm1out_000001.nc')
# ds = nc.Dataset('C:/Users/mschne28/Documents/cm1r21.1/run/cm1out_000001.nc')
zh = ds.variables['zh'][:].data
u0 = ds.variables['uinterp'][:].data[0,:,10,10]
v0 = ds.variables['vinterp'][:].data[0,:,10,10]
# u0 = ds.variables['u0'][:].data[0,:,0,0]
# v0 = ds.variables['v0'][:].data[0,:,0,0]
umove = ds.variables['umove'][:].data[0]
vmove = ds.variables['vmove'][:].data[0]
ds.close()

u0 = u0 + umove
v0 = v0 + vmove

i1 = np.argmin(abs(zh-0.5))
i2 = np.argmin(abs(zh-5.5))
i3 = np.argmin(abs(zh-6))

D = 7.5
u_06 = np.mean(u0[0:81])
v_06 = np.mean(v0[0:81])
u_shr = np.mean(u0[77:81]) - np.mean(u0[0:25])
v_shr = np.mean(v0[77:81]) - np.mean(v0[0:25])
V_shr = np.sqrt(u_shr**2 + v_shr**2)
u_rm = u_06 + D*(v_shr/V_shr)
v_rm = v_06 - D*(u_shr/V_shr)


phi = np.linspace(0,360,361)

fig,ax = plt.subplots(1, 1, figsize=(5,5), subplot_kw=dict(box_aspect=1), layout='constrained')

ax.axhline(0, color='grey', linewidth=1)
ax.axvline(0, color='grey', linewidth=1)
ax.plot(10*np.cos(phi*np.pi/180), 10*np.sin(phi*np.pi/180), 'lightgrey', linewidth=1)
ax.plot(20*np.cos(phi*np.pi/180), 20*np.sin(phi*np.pi/180), 'lightgrey', linewidth=1)
ax.plot(30*np.cos(phi*np.pi/180), 30*np.sin(phi*np.pi/180), 'lightgrey', linewidth=1)
ax.plot(40*np.cos(phi*np.pi/180), 40*np.sin(phi*np.pi/180), 'lightgrey', linewidth=1)
ax.plot(50*np.cos(phi*np.pi/180), 50*np.sin(phi*np.pi/180), 'lightgrey', linewidth=1)
ax.plot(u0, v0, 'k', linewidth=2)
ax.plot([u0[0],u_06], [v0[0],v_06], 'k', linewidth=3)
ax.plot([u0[0],u_shr], [v0[0],v_shr], 'b', linewidth=3)
ax.plot([u0[0],u_rm], [v0[0],v_rm], 'r', linewidth=3)
ax.set_title('Base state hodograph', fontsize=10)
ax.set_xlim([-40,40])
ax.set_ylim([-40,40])

plt.show()




#%% more cold blob mods/debugging

fp = 'C:/Users/mschne28/Documents/cm1out/warmblob_test/'
ds = nc.Dataset(fp+'cm1out_000001.nc')
xh = ds.variables['xh'][:].data*1000
yh = ds.variables['yh'][:].data*1000
zh = ds.variables['zh'][:].data*1000
ds.close()

# xh = np.linspace(-60000,60000,121)
# yh = np.linspace(-60000,60000,121)

s0 = -0.020 # max nudging
t0 = 0
t1 = 900 #end rampup
t2 = 901 #start rampdown
t3 = 1800 #end rampdown

thp = (t3+t2-t1)*s0/2





t = np.linspace(0,1800,1801)
s = np.zeros(shape=t.shape, dtype=float)
# th = np.zeros(shape=t.shape, dtype=float)

s[(t<t1)] = s0*t[(t<t1)]/t1
s[(t>=t1)&(t<t2)] = s0
s[(t>=t2)&(t<t3)] = s0*(t3-t[(t>=t2)&(t<t3)])/(t3-t2)

# th = [np.sum(s[:n]) for n in range(len(s))]


coldmod = np.zeros(shape=(len(t),), dtype=float)
bigR = np.zeros(shape=(len(yh),len(xh)), dtype=float)
bigZ = np.zeros(shape=(len(zh),), dtype=float)

coldxc = 0
coldyc = 0
coldzc = 0
coldrad = 10000
coldrady = 40000
colddepth = 3000
coldamp = s0


for l in range(len(t)):
    if t[l] < t1:
        coldmod[l] = (t[l]-t0)/(t1-t0)
    elif (t[l] >= t1) & (t[l] < t2):
        coldmod[l] = 1
    elif (t[l] >= t2) & (t[l] < t3):
        coldmod[l] = (t3-t[l])/(t3-t2)
    else:
        coldmod[l] = 0

for j in range(len(yh)):
    for i in range(len(xh)):
        if coldrady <= 0:
            radius = ((xh[i]-coldxc)**2 + (yh[j]-coldyc)**2)**0.5
            if radius < coldrad:
                bigR[j,i] = 1 - (radius**2/coldrad**2)
            else:
                bigR[j,i] = 0
        else:
            xrad = np.abs(xh[i] - coldxc)
            yrad = np.abs(yh[j] - coldyc)
            crad = coldrad + coldrady
            if xrad < coldrad:
                if yrad < coldrady:
                    radius = np.abs(xh[i]-coldxc)
                    bigR[j,i] = 1 - (radius**2 / coldrad**2)
                elif ((yrad >= coldrady) & (yrad < crad)):
                    radius = ((xh[i]-coldxc)**2 + (np.abs(yh[j]-coldyc)-coldrady)**2)**0.5
                    bigR[j,i] = 1 - (radius**2 / coldrad**2)
                else:
                    bigR[j,i] = 0
            else:
                bigR[j,i] = 0

for k in range(len(zh)):
    radius = np.abs(zh[k] - coldzc)
    if radius < colddepth:
        bigZ[k] = 1 - (radius**2/colddepth**2)
    else:
        bigZ[k] = 0



# fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
# plot_contourf(xh, yh, bigR, 'dbz', ax, levels=np.linspace(0,1,21),
#               datalims=[0,1], cmap='Reds')
# ax.set_title(f"bigR")
# plt.show()



ll = 150
xl = [-ll,ll]
yl = [-ll,ll]

di = np.linspace(0,1800,5)
iz = np.argmin(abs(zh - 0))
for i in di:
    # S = bigR*s[int(i)]
    # th = bigR*np.sum(s[:int(i)])
    S = bigR*bigZ[iz]*coldmod[int(i)]*coldamp
    th = bigR*bigZ[iz]*np.sum(coldmod[:int(i)])*coldamp
    
    x = 0.001*xh
    y = 0.001*yh
    z = 0.001*zh
    
    
    fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
    # plot_contourf(x, y, S, 'dbz', ax, levels=np.linspace(0,0.02,21),
    #               datalims=[0,0.02], xlims=xl, ylims=yl, cmap='Reds')
    plot_contourf(x, y, S, 'dbz', ax, levels=np.linspace(-0.02,0,21),
                  datalims=[-0.02,0], xlims=xl, ylims=yl, cmap='Blues_r')
    ax.set_title(f"S, t = {i} s")
    plt.show()
    
    
    # fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
    # # plot_contourf(x, y, th, 'th', ax, levels=np.linspace(0,18,19),
    # #               datalims=[0,18], xlims=xl, ylims=yl, cmap='Reds')
    # plot_contourf(x, y, th, 'th', ax, levels=np.linspace(-18,0,19),
    #               datalims=[-18,0], xlims=xl, ylims=yl, cmap='Blues_r')
    # ax.set_title(f"theta, t = {i} s")
    # # plt.show()
    


bigZ2 = np.tile(bigZ, [len(yh),len(xh),1]).transpose()
bigR2 = np.tile(bigR, [len(zh),1,1])

sk = bigZ2[:,:,0] * bigR2[:,:,149] * coldmod[1200] * coldamp


fig,ax = plt.subplots(2, 2, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')

plot_contourf(xh, yh, bigR, 'dbz', ax[0,0], levels=np.linspace(0,1,21),
              datalims=[0,1], cmap='Reds')
ax[0,0].set_title(f"bigR")

plot_contourf(yh, zh, np.tile(bigZ, [len(yh),1]).transpose(), 'dbz', ax[0,1],
              levels=np.linspace(0,1,21), datalims=[0,1], cmap='Reds')
# ax[0,1].plot(zh, bigZ, 'k', linewidth=2)
# ax[0,1].set_ylim([0,1])
ax[0,1].set_title(f"bigZ")

ax[1,0].plot(t, coldmod, 'k', linewidth=2)
ax[1,0].set_ylim([0,1])
ax[1,0].set_title(f"coldmod")

plot_contourf(yh, zh, sk, 'dbz', ax[1,1], levels=np.linspace(-0.02,0,21),
                 datalims=[-0.02,0], cmap='Blues_r')
ax[1,1].set_title(f"S")

plt.show()



fig,ax = plt.subplots(1, 1, figsize=(8,4))
ax.plot(t, coldmod, 'k', linewidth=2)
plt.show

fig,ax = plt.subplots(1, 1, figsize=(8,4))
ax.plot(np.min(bigZ, axis=(1,2)), zh, 'k', linewidth=2)
plt.show


#%% recalculating dbz components to check why p3 is weird

dblmom = False

# ds = nc.Dataset('C:/Users/mschne28/Documents/cm1r21.1/run/cm1out_000007.nc')
ds = nc.Dataset('C:/Users/mschne28/Documents/cm1out/wk_p3_500m/cm1out_000022.nc')

rho_l = 1000
rho_i = 900
rho = 1.1

time = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['zh'][:].data

dbz = ds.variables['dbz'][:].data[0,:,:,:]
cref = ds.variables['cref'][:].data[0,:,:]
qc = ds.variables['qc'][:].data[0,:,:,:] #kg/kg
qr = ds.variables['qr'][:].data[0,:,:,:]
qi = ds.variables['qi'][:].data[0,:,:,:]
if dblmom:
    qi2 = ds.variables['qi2'][:].data[0,:,:,:]
    qir = ds.variables['ri'][:].data[0,:,:,:]
    qir2 = ds.variables['ri2'][:].data[0,:,:,:]
    qnc = ds.variables['nc'][:].data[0,:,:,:]
    qnr = ds.variables['nr'][:].data[0,:,:,:]
    qni = ds.variables['ni'][:].data[0,:,:,:]
    qni2 = ds.variables['ni2'][:].data[0,:,:,:]
    qbi = ds.variables['bi'][:].data[0,:,:,:] #1/m3/kg rime ice volume
    qbi2 = ds.variables['bi2'][:].data[0,:,:,:]
else:
    qir = ds.variables['qir'][:].data[0,:,:,:] #rime ice mass mixing ratio kg/kg
    qnc = ds.variables['qnc'][:].data[0,:,:,:] #1/kg number conc
    qnr = ds.variables['qnr'][:].data[0,:,:,:]
    qni = ds.variables['qni'][:].data[0,:,:,:]
    qbi = ds.variables['qib'][:].data[0,:,:,:] #1/m3/kg rime ice volume
    qzi = ds.variables['qzi'][:].data[0,:,:,:] #1/m6/kg ice reflectivity

ds.close()


zr = np.where(qnr>0, 20*1e18*(6/(np.pi*rho_l))**2 * (rho*qr)**2 / qnr, 0)


g0 = (6*5*4)/(3*2*1)
g05 = (6.5*5.5*4.5)/(3.5*2.5*1.5)
g20 = (26*25*24)/(23*22*21)

if dblmom:
    zi1 = np.where(qni>0, 0.224*1e18*g05*(6/(np.pi*rho_i))**2 * (rho*qi)**2/qni, 0)
    zi2 = np.where(qni2>0, 0.224*1e18*g05*(6/(np.pi*rho_i))**2 * (rho*qi2)**2/qni2, 0)
    
    z01 = np.where(qni>0, 0.224*1e18*g0*(6/(np.pi*rho_i))**2 * (rho*qi)**2/qni, 0)
    z201 = np.where(qni>0, 0.224*1e18*g20*(6/(np.pi*rho_i))**2 * (rho*qi)**2/qni, 0)
    z02 = np.where(qni2>0, 0.224*1e18*g0*(6/(np.pi*rho_i))**2 * (rho*qi2)**2/qni2, 0)
    z202 = np.where(qni2>0, 0.224*1e18*g20*(6/(np.pi*rho_i))**2 * (rho*qi2)**2/qni2, 0)
    
    zi1 = np.minimum(zi1, z01)
    zi1 = np.maximum(zi1, z201)
    zi2 = np.minimum(zi2, z02)
    zi2 = np.maximum(zi2, z202)
    
    zi = zi1 + zi2
    
    Zi1 = np.where(zi1>0, 10*np.log10(zi1), 0)
    Zi2 = np.where(zi2>0, 10*np.log10(zi2), 0)
    
    Zic1 = np.max(Zi1, axis=0)
    Zic2 = np.max(Zi2, axis=0)
else:
    zi = np.where(qni>0, 0.224*1e18*g05*(6/(np.pi*rho_i))**2 * (rho*qi)**2/qni, 0)
    z0 = np.where(qni>0, 0.224*1e18*g0*(6/(np.pi*rho_i))**2 * (rho*qi)**2/qni, 0)
    z20 = np.where(qni>0, 0.224*1e18*g20*(6/(np.pi*rho_i))**2 * (rho*qi)**2/qni, 0)
    
    zi = np.minimum(zi, z0)
    zi = np.maximum(zi, z20)


Zi = np.where(zi>0, 10*np.log10(zi), 0)
Zr = np.where(zr>0, 10*np.log10(zr), 0)
Ze = np.where(zr+zi>0, 10*np.log10(zr+zi), 0)

Zrc = np.max(Zr, axis=0)
Zic = np.max(Zi, axis=0)
Zec = np.max(Ze, axis=0)

fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, np.ma.masked_array(Zrc, Zrc<0.1), 'dbz', ax, levels=np.linspace(0,70,15), datalims=[0,70], cmap='HomeyerRainbow')
ax.set_title('rain comp dbz (calculated)')


fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, np.ma.masked_array(Zic, Zic<0.1), 'dbz', ax, levels=np.linspace(0,70,15), datalims=[0,70], cmap='HomeyerRainbow')
ax.set_title('ice comp dbz (calculated)')

if dblmom:
    fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
    plot_contourf(xh, yh, np.ma.masked_array(Zic1, Zic1<0.1), 'dbz', ax, levels=np.linspace(0,70,15), datalims=[0,70], cmap='HomeyerRainbow')
    ax.set_title('ice 1 comp dbz (calculated)')
    
    fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
    plot_contourf(xh, yh, np.ma.masked_array(Zic2, Zic2<0.1), 'dbz', ax, levels=np.linspace(0,70,15), datalims=[0,70], cmap='HomeyerRainbow')
    ax.set_title('ice 2 comp dbz (calculated)')


fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, np.ma.masked_array(Zec, Zec<0.1), 'dbz', ax, levels=np.linspace(0,70,15), datalims=[0,70], cmap='HomeyerRainbow')
ax.set_title('total comp dbz (calculated)')

fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
plot_contourf(xh, yh, np.ma.masked_array(cref, cref<0.1), 'dbz', ax, levels=np.linspace(0,70,15), datalims=[0,70], cmap='HomeyerRainbow')
ax.set_title('cref (model)')




#%% resave all the cm1out_stats variables for the noslip sim because I fucked up

# Reprocess reflectivity if needed
if False:
    ds = nc.Dataset('C:/Users/mschne28/Documents/cm1out/noslip_wk_250m/cm1out_stats.nc')
    ds2 = nc.Dataset('C:/Users/mschne28/Documents/cm1out/noslip_wk_250m/cm1out_stats2.nc', 'r+', clobber=True)
    
    xh = ds.variables['xh'][:].data
    yh = ds.variables['yh'][:].data
    zh = ds.variables['zh'][:].data
    newtime = np.linspace(0,25200,421)
    
    ds2.createDimension('xh', size=len(xh))
    ds2.createDimension('yh', size=len(yh))
    ds2.createDimension('zh', size=len(zh))
    ds2.createDimension('time', size=len(newtime))
    
    cx = ds2.createVariable('xh', 'f4', ('xh'))
    cx.units = 'm'; cx.long_name = 'west-east location'; cx.axis = 'X'; cx[:] = xh[:]
    cy = ds2.createVariable('yh', 'f4', ('yh'))
    cy.units = 'm'; cy.long_name = 'south-north location'; cy.axis = 'Y'; cy[:] = yh[:]
    cz = ds2.createVariable('zh', 'f4', ('zh'))
    cz.units = 'm'; cz.long_name = 'height'; cz.axis = 'Z'; cz[:] = zh[:]
    ct = ds2.createVariable('time', 'f4', ('time'))
    ct.units = 's'; ct.long_name = 'time'; ct.axis='T'; ct[:] = newtime[:]
    cmt = ds2.createVariable('mtime', 'f4', ('time'))
    cmt.units = 's'; cmt.long_name = 'model time (seconds since beginning of s)'; cmt[:] = newtime[:]
    
    dvars = list(ds.variables.keys())[5:]
    for v in dvars:
        print(f"Variable: {v}")
        data = ds.variables[v][:].data[:421]
        tmp = ds2.createVariable(v, 'f4', ('time'))
        tmp.units = ds.variables[v].units
        tmp.long_name = ds.variables[v].long_name
        tmp[:] = data[:]
    
    ds2.close()
    ds.close()





