# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 17:00:12 2026

@author: mschne28
"""

from CM1utils import *

#%% Load data?

fp = 'C:/Users/mschne28/Documents/cm1out/brooks/era5-1_125m_test5/'

ds = nc.Dataset(fp+f"cm1out_000005.nc")

time = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['zh'][:].data

# get z indices
ix1 = np.argmin(abs(xh+50))
ix2 = np.argmin(abs(xh-50))
iy1 = np.argmin(abs(yh+50))
iy2 = np.argmin(abs(yh-50))
ix = slice(ix1,ix2+1)
iy = slice(iy1,iy2+1)

xf = xh[ix]
yf = yh[iy]

iz2 = np.argmin(abs(zh-2))
iz80 = np.argmin(abs(zh-0.07))
iz05 = np.argmin(abs(zh-0.5))
iz1 = np.argmin(abs(zh-1))

# winterp = ds.variables['winterp'][:].data[0,:iz2+2,:,:]
# uinterp = ds.variables['uinterp'][:].data[0,:iz2+2,:,:] + ds.variables['umove'][:].data[0]
# vinterp = ds.variables['vinterp'][:].data[0,:iz2+2,:,:] + ds.variables['vmove'][:].data[0]
# u80 = np.mean(uinterp[iz80:iz80+2,:,:], axis=0)
# v80 = np.mean(vinterp[iz80:iz80+2,:,:], axis=0)
# V_80m = np.sqrt(u80**2 + v80**2)
# V_2km = np.sqrt(uinterp[iz2,:,:]**2 + vinterp[iz2,:,:]**2)
# w_2km = np.mean(winterp[iz05:iz2+1,:,:], axis=0)

dbz = np.mean(ds.variables['dbz'][:].data[0,iz80:iz80+2,:,:], axis=0)

# 80-m wind criteria
u80m = np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:] + ds.variables['umove'][:].data[0], axis=0)
v80m = np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:] + ds.variables['vmove'][:].data[0], axis=0)
V80m = np.sqrt(u80m**2 + v80m**2)

# RIJ criteria
u2km = ds.variables['uinterp'][:].data[0,iz2,:,:] + ds.variables['umove'][:].data[0]
v2km = ds.variables['vinterp'][:].data[0,iz2,:,:] + ds.variables['vmove'][:].data[0]
V2km = np.sqrt(u2km**2 + v2km**2)
w2km = np.mean(ds.variables['winterp'][:].data[0,0:iz2+1,:,:], axis=0)

# MV criteria? what height? not specified in Killian thesis - check Lasher-Trapp et al. 2023 *** i guess it's at 80 m?
zvort = np.mean(ds.variables['zvort'][:].data[0,iz80:iz80+2,:,:], axis=0)
D1 = np.gradient(u80m, xh*1000, axis=1) - np.gradient(v80m, yh*1000, axis=0)
D2 = np.gradient(v80m, xh*1000, axis=1) + np.gradient(u80m, yh*1000, axis=0)
OW80m = zvort**2 - D1**2 - D2**2

# DB criteria?
w1km = ds.variables['winterp'][:].data[0,iz1,:,:]
w_dn_max = np.max(ds.variables['winterp'][:].data[0,:,:,:], axis=0)

ds.close()


#%


# Conditions from Killian thesis - adapted from Lasher-Trapp et al. 2023
Vsub_thres = 20 # sub-svr wind threshold
V_thres = 25.7 #wind speed threshold
Vsig_thres = 33.4 # sig-svr wind threshold
w_thres_rij = -2 #downdraft threshold
w_thres_db = -5
ow_thres_mv = 1e-4

is_rij = np.zeros(shape=(len(yh),len(xh)), dtype=int)
V2_flag = np.zeros(shape=(len(yh),len(xh)), dtype=int) # Killian thesis - use 2 km wspd for RIJ ID
w2_flag = np.zeros(shape=(len(yh),len(xh)), dtype=int) # 0.5-2 km or 0-2 km? unclear from Killian thesis
sub_flag = np.zeros(shape=(len(yh),len(xh)), dtype=int) # flag sub-svr 80-m wind
svr_flag = np.zeros(shape=(len(yh),len(xh)), dtype=int) # flag svr 80-m wind
sig_flag = np.zeros(shape=(len(yh),len(xh)), dtype=int) # flag sig-svr 80-m wind
is_mv = np.zeros(shape=(len(yh),len(xh)), dtype=int)
ow_flag = np.zeros(shape=(len(yh),len(xh)), dtype=int)
is_db = np.zeros(shape=(len(yh),len(xh)), dtype=int)
w1_flag = np.zeros(shape=(len(yh),len(xh)), dtype=int)
wmax_flag = np.zeros(shape=(len(yh),len(xh)), dtype=int)
is_mv_rij = np.zeros(shape=(len(yh),len(xh)), dtype=int)
is_mv_db = np.zeros(shape=(len(yh),len(xh)), dtype=int)

# Svr/sig svr wind ID
sub_flag[(V80m >= Vsub_thres)] = 1
svr_flag[(V80m >= V_thres)] = 1
sig_flag[(V80m >= Vsig_thres)] = 1

# RIJ ID
V2_flag[(V2km >= V_thres)] = 1
w2_flag[(w2km <= w_thres_rij)] = 1

# MV ID
ow_flag[(OW80m >= ow_thres_mv)] = 1

# DB ID
w1_flag[(w1km <= w_thres_db)] = 1
wmax_flag[(w_dn_max < 0)] = 1


# Find criteria in a 5 km x 5 km box around each point in the 125-m subdomain
for j in range(len(yf)):
    for i in range(len(xf)):
        idx1 = np.argmin(abs(xh-(xf[i]-5)))
        idx2 = np.argmin(abs(xh-(xf[i]+5)))
        idy1 = np.argmin(abs(yh-(yf[j]-5)))
        idy2 = np.argmin(abs(yh-(yf[j]+5)))
        idx = slice(idx1,idx2+1)
        idy = slice(idy1,idy2+1)
        ixc = np.argmin(abs(xh-xf[i]))
        iyc = np.argmin(abs(yh-yf[j]))
        
        # RIJ criteria
        if (np.max(svr_flag[idy,idx]) > 0) & (np.max(V2_flag[idy,idx]) > 0) & (np.max(w2_flag[idy,idx]) > 0):
            is_rij[iyc,ixc] = 1
        
        # MV criteria
        if (np.max(svr_flag[idy,idx]) > 0) & (np.max(ow_flag[idy,idx]) > 0):
            is_mv[iyc,ixc] = 1
        
        # DB criteria
        if (np.max(svr_flag[idy,idx]) > 0) & (np.max(w1_flag[idy,idx]) > 0) & (np.max(wmax_flag[idy,idx]) > 0):
            is_db[iyc,ixc] = 1
        
        # MV+RIJ criteria
        if (is_mv[iyc,ixc]) & (is_rij[iyc,ixc]):
            is_mv_rij[iyc,ixc] = 1
        
        # MV+DB criteria
        if (is_mv[iyc,ixc]) & (is_db[iyc,ixc]):
            is_mv_db[iyc,ixc] = 1


# is_rij[(svr_flag) & (V2_flag) & (w2_flag)] = 1


if True:
    dbfile = open(fp+f"wind_mechanisms_{time/60:.0f}min.pkl", 'wb')
    data = {'is_rij':is_rij, 'is_mv':is_mv, 'is_db':is_db, 'is_mv_rij':is_mv_rij, 'is_mv_db':is_mv_db}
    pickle.dump(data, dbfile)
    dbfile.close()


if False:
    fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1))
    plot_cfill(xh, yh, np.ma.masked_array(sub_flag, sub_flag<1), 'w', ax, datalims=[0,1], cmap='Greys')
    plot_cfill(xh, yh, np.ma.masked_array(svr_flag, svr_flag<1), 'w', ax, datalims=[0,1], cmap='managua', cbar=False)
    plot_cfill(xh, yh, np.ma.masked_array(sig_flag, sig_flag<1), 'w', ax, datalims=[0,1], cmap='BuOrR14', cbar=False)
    ax.set_xlim([-50,50])
    ax.set_ylim([-50,50])
    ax.set_title(f"Sub-severe, severe, and sig-severe 80-m wind (t = {time/60:.0f} min)", fontsize=12)
    plt.show()
    
    
    
    fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1))
    plot_cfill(xh, yh, np.ma.masked_array(w_flag, w_flag<1), 'w', ax, datalims=[0,1], cmap='Bu10')
    plot_cfill(xh, yh, np.ma.masked_array(V_flag, V_flag<1), 'wspd', ax, datalims=[0,1], cmap='BuOrR14', cbar=False)
    ax.set_xlim([-50,50])
    ax.set_ylim([-50,50])
    ax.set_title(f"RIJ criteria - w < -2 m/s, V > 20 m/s (t = {time/60:.0f} min)", fontsize=12)
    plt.show()
    
#%%
    
# fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1))
# plot_cfill(xh, yh, np.ma.masked_array(is_rij, is_rij<1), 'w', ax, datalims=[0,1], cmap='Greys')
# ax.contour(xh, yh, w2_flag, levels=[0.1], colors='deepskyblue', linewidths=1)
# ax.contour(xh, yh, V2_flag, levels=[0.1], colors='gold', linewidths=1.5)
# ax.contour(xh, yh, svr_flag, levels=[0.1], colors='r', linewidths=1.5)
# ax.set_xlim([-50,50])
# ax.set_ylim([-50,50])
# ax.set_title(f"RIJ criteria (t = {time/60:.0f} min)", fontsize=12)
# l1, = ax.plot([149,150], [149,150], 'deepskyblue', linewidth=1)
# l2, = ax.plot([149,150], [149,150], 'gold', linewidth=1.5)
# l3, = ax.plot([149,150], [149,150], 'r', linewidth=1.5)
# ax.legend(handles=[l1,l2,l3], labels=['w(0-2 km) < -2 m/s','V(2 km) > 25.7 m/s', 'V(80 m) > 25.7 m/s'], loc='lower left', fontsize=10)
# plt.show()



rij_mask = (is_rij<1) | (is_mv_rij)
mv_mask = (is_mv<1) | (is_mv_rij) | (is_mv_db)
db_mask = (is_db<1) | (is_mv_db)
mv_rij_mask = (is_mv_rij<1)
mv_db_mask = (is_mv_db<1)

fig,ax = plt.subplots(1, 1, figsize=(9,6), subplot_kw=dict(box_aspect=1))

plot_contourf(xh, yh, np.ma.masked_array(dbz, dbz<10), 'dbz', ax, levels=np.linspace(0,70,15), datalims=[0,70])

# plot_cfill(xh, yh, np.ma.masked_array(is_rij, is_rij<1), 'w', ax, datalims=[0,1], cmap='PiYG', cbar=False)
# # plot_cfill(xh, yh, np.ma.masked_array(is_mv, is_mv<1), 'w', ax, datalims=[0,1], cmap='bwr', cbar=False)
# # plot_cfill(xh, yh, np.ma.masked_array(is_db, is_db<1), 'w', ax, datalims=[0,1], cmap='Bu10', cbar=False)
# # plot_cfill(xh, yh, np.ma.masked_array(is_mv_rij, is_mv_rij<1), 'w', ax, datalims=[0,1], cmap='spring_r', cbar=False)
# # plot_cfill(xh, yh, np.ma.masked_array(is_mv_db, is_mv_db<1), 'w', ax, datalims=[0,1], cmap='managua_r', cbar=False)
# # ax.contour(xh, yh, is_rij, levels=[0.1], colors='green', linewidths=1.5)
# ax.contour(xh, yh, is_mv, levels=[0.9], colors='r', linewidths=1.5)
# ax.contour(xh, yh, is_db, levels=[0.9], colors='dodgerblue', linewidths=1.5)
# # ax.contour(xh, yh, np.ma.masked_array(is_mv, (is_mv_rij) | (is_mv_db)), levels=[0.9], colors='r', linewidths=1.5)
# # ax.contour(xh, yh, np.ma.masked_array(is_db, is_mv_db), levels=[0.9], colors='dodgerblue', linewidths=1.5)
# ax.contour(xh, yh, is_mv_rij, levels=[0.9], colors='violet', linewidths=2)
# ax.contour(xh, yh, is_mv_db, levels=[0.9], colors='gold', linewidths=1.5)

plot_cfill(xh, yh, np.ma.masked_array(is_rij, rij_mask), 'w', ax, datalims=[0,1], cmap='PiYG', cbar=False)
plot_cfill(xh, yh, np.ma.masked_array(is_mv, mv_mask), 'w', ax, datalims=[0,1], cmap='bwr', cbar=False)
plot_cfill(xh, yh, np.ma.masked_array(is_db, db_mask), 'w', ax, datalims=[0,1], cmap='Bu10', cbar=False)
plot_cfill(xh, yh, np.ma.masked_array(is_mv_rij, mv_rij_mask), 'w', ax, datalims=[0,1], cmap='vanimo_r', cbar=False)
plot_cfill(xh, yh, np.ma.masked_array(is_mv_db, mv_db_mask), 'w', ax, datalims=[0,1], cmap='managua_r', cbar=False)

ax.contour(xh, yh, V80m, levels=[25], colors='k', linewidths=[1.5])

ax.set_xlim([-50,50])
ax.set_ylim([-50,50])
ax.set_title(f"Severe wind mechanisms (t = {time/60:.0f} min)", fontsize=12)
l1, = ax.plot([149,150], [149,150], 'green', linewidth=1.5)
l2, = ax.plot([149,150], [149,150], 'r', linewidth=1.5)
l3, = ax.plot([149,150], [149,150], 'dodgerblue', linewidth=1.5)
l4, = ax.plot([149,150], [149,150], 'violet', linewidth=1.5)
l5, = ax.plot([149,150], [149,150], 'gold', linewidth=1.5)
l6, = ax.plot([149,150], [149,150], 'k', linewidth=1.5)
ax.legend(handles=[l1,l2,l3,l4,l5,l6], labels=['RIJ','MV','DB','MV+RIJ','MV+DB','SVR'], loc='lower left', fontsize=10)
plt.show()


# fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1))
# plot_contourf(xh, yh, np.ma.masked_array(w2km, w2km > w_thres_rij), 'w', ax, levels=np.linspace(-6,-2,11), datalims=[-6,0], cmap='Bu10_r')
# ax.contour(xh, yh, V2_flag, levels=[0.1], colors='r', linewidths=1)
# ax.set_xlim([-50,50])
# ax.set_ylim([-50,50])
# ax.set_title(f"w < -2 m/s, V > 20 m/s at 2 km (t = {time/60:.0f} min)", fontsize=12)
# plt.show()




#%% Add wind mechanisms to translated swaths

fp = 'C:/Users/mschne28/Documents/cm1out/brooks/era5-1_125m_test5/'



ds = nc.Dataset(fp+'cm1out_000005.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
umove = ds.variables['umove'][:].data[0]
vmove = ds.variables['vmove'][:].data[0]
sws1 = ds.variables['sws2'][:].data[0,:,:]
shs1 = ds.variables['shs2'][:].data[0,:,:]
hail1 = ds.variables['hail2'][:].data[0,:,:]
dbz1 = ds.variables['dbz'][:].data[0,0,:,:]
wsp1 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove, axis=0))**2 + 
               (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove, axis=0))**2)
ds.close()

ds = nc.Dataset(fp+'cm1out_000009.nc')
sws2 = ds.variables['sws2'][:].data[0,:,:]
shs2 = ds.variables['shs2'][:].data[0,:,:]
hail2 = ds.variables['hail2'][:].data[0,:,:]
dbz2 = ds.variables['dbz'][:].data[0,0,:,:]
wsp2 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove, axis=0))**2 + 
               (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove, axis=0))**2)
ds.close()

ds = nc.Dataset(fp+'cm1out_000013.nc')
sws3 = ds.variables['sws2'][:].data[0,:,:]
shs3 = ds.variables['shs2'][:].data[0,:,:]
hail3 = ds.variables['hail2'][:].data[0,:,:]
dbz3 = ds.variables['dbz'][:].data[0,0,:,:]
wsp3 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove, axis=0))**2 + 
               (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove, axis=0))**2)
ds.close()

ds = nc.Dataset(fp+'cm1out_000017.nc')
sws4 = ds.variables['sws2'][:].data[0,:,:]
shs4 = ds.variables['shs2'][:].data[0,:,:]
hail4 = ds.variables['hail2'][:].data[0,:,:]
dbz4 = ds.variables['dbz'][:].data[0,0,:,:]
wsp4 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove, axis=0))**2 + 
               (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove, axis=0))**2)
ds.close()

ds = nc.Dataset(fp+'cm1out_000021.nc')
sws5 = ds.variables['sws2'][:].data[0,:,:]
shs5 = ds.variables['shs2'][:].data[0,:,:]
hail5 = ds.variables['hail2'][:].data[0,:,:]
dbz5 = ds.variables['dbz'][:].data[0,0,:,:]
wsp5 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove, axis=0))**2 + 
               (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove, axis=0))**2)
ds.close()

ds = nc.Dataset(fp+'cm1out_000025.nc')
sws6 = ds.variables['sws2'][:].data[0,:,:]
shs6 = ds.variables['shs2'][:].data[0,:,:]
hail6 = ds.variables['hail2'][:].data[0,:,:]
dbz6 = ds.variables['dbz'][:].data[0,0,:,:]
wsp6 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove, axis=0))**2 + 
               (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove, axis=0))**2)
ds.close()

ds = nc.Dataset(fp+'cm1out_000029.nc')
sws7 = ds.variables['sws2'][:].data[0,:,:]
shs7 = ds.variables['shs2'][:].data[0,:,:]
hail7 = ds.variables['hail2'][:].data[0,:,:]
dbz7 = ds.variables['dbz'][:].data[0,0,:,:]
wsp7 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove, axis=0))**2 + 
               (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove, axis=0))**2)
ds.close()

ds = nc.Dataset(fp+'cm1out_000033.nc')
sws8 = ds.variables['sws2'][:].data[0,:,:]
shs8 = ds.variables['shs2'][:].data[0,:,:]
hail8 = ds.variables['hail2'][:].data[0,:,:]
dbz8 = ds.variables['dbz'][:].data[0,0,:,:]
wsp8 = np.sqrt((np.mean(ds.variables['uinterp'][:].data[0,iz80:iz80+2,:,:]+umove, axis=0))**2 + 
               (np.mean(ds.variables['vinterp'][:].data[0,iz80:iz80+2,:,:]+vmove, axis=0))**2)
ds.close()



dbfile = open(fp+'wind_mechanisms_60min.pkl', 'rb')
crit1 = pickle.load(dbfile)
dbfile.close()

dbfile = open(fp+'wind_mechanisms_120min.pkl', 'rb')
crit2 = pickle.load(dbfile)
dbfile.close()

dbfile = open(fp+'wind_mechanisms_180min.pkl', 'rb')
crit3 = pickle.load(dbfile)
dbfile.close()

dbfile = open(fp+'wind_mechanisms_240min.pkl', 'rb')
crit4 = pickle.load(dbfile)
dbfile.close()

dbfile = open(fp+'wind_mechanisms_300min.pkl', 'rb')
crit5 = pickle.load(dbfile)
dbfile.close()

dbfile = open(fp+'wind_mechanisms_360min.pkl', 'rb')
crit6 = pickle.load(dbfile)
dbfile.close()

dbfile = open(fp+'wind_mechanisms_420min.pkl', 'rb')
crit7 = pickle.load(dbfile)
dbfile.close()

dbfile = open(fp+'wind_mechanisms_480min.pkl', 'rb')
crit8 = pickle.load(dbfile)
dbfile.close()

crit = {'1':crit1, '2':crit2, '3':crit3, '4':crit4, '5':crit5, '6':crit6, '7':crit7, '8':crit8}


x_added = umove*3600/1000
y_added = vmove*3600/1000

xh1 = xh + 100
xh2 = xh1 + x_added
xh3 = xh1 + 2*x_added
xh4 = xh1 + 3*x_added
xh5 = xh1 + 4*x_added
xh6 = xh1 + 5*x_added
xh7 = xh1 + 6*x_added
xh8 = xh1 + 7*x_added

yh1 = yh + 100
yh2 = yh1 + y_added
yh3 = yh1 + 2*y_added
yh4 = yh1 + 3*y_added
yh5 = yh1 + 4*y_added
yh6 = yh1 + 5*y_added
yh7 = yh1 + 6*y_added
yh8 = yh1 + 7*y_added

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
    xt = [77, 120, 165, 215, 265, 320, 380, 435]
    yt = [73, 78, 80, 83, 89, 96, 102, 104]
elif 'era5-1_125m_test2' in fp:
    xt = [77, 120, 165, 215, 265, 325, 378, 425]
    yt = [70, 74, 77, 82, 87, 93, 94, 96]
elif 'era5-1_125m_test3' in fp:
    xt = [77, 116, 160, 210, 255, 306, 355, 405]
    yt = [71, 75, 79, 79, 80, 81, 82, 86]
elif 'era5-1_125m_test4' in fp:
    xt = [75, 115, 157, 201, 247, 298, 348, 397]
    yt = [74, 83, 85, 85, 85, 85, 85, 85]
elif 'era5-1_125m_test5' in fp:
    xt = [78, 120, 160, 202, 250, 295, 335, 380]
    yt = [65, 73, 73, 75, 75, 78, 79, 79]
elif 'era5-1_125m_test6' in fp:
    xt = [80, 122, 163, 205, 250, 292, 338, 380]
    yt = [60, 67, 68, 70, 71, 73, 73, 73]
elif 'era5-1_125m_test7' in fp:
    xt = [78, 120, 165, 208, 255, 300, 342, 383]
    yt = [65, 68, 72, 73, 75, 75, 73, 72]
    
    
xl = [50,500]
yl = [50,200]


if 'era5-1' in fp:
    sounding_str = 'ERA5 profile (50.75,-114.0)'
elif 'era5-2' in fp:
    sounding_str = 'ERA5 profile (51.0,-114.25)'
elif 'era5-3' in fp:
    sounding_str = 'ERA5 profile (51.25,-113.25)'
elif 'hrdps' in fp:
    sounding_str = 'HRDPS profile (51.166,-113.135)'



figsave = False

levs = np.linspace(0,30,31)
cm = 'Blues'


fig,ax = plt.subplots(1, 1, figsize=(8.5,2.75), subplot_kw=dict(aspect=1), layout='constrained')

c = ax.contourf(xh1, yh1, np.ma.masked_array(wsp1, dbz1<20), levels=levs, vmin=0, vmax=30, cmap=cm)
ax.contourf(xh2, yh2, np.ma.masked_array(wsp2, dbz2<20), levels=levs, vmin=0, vmax=30, cmap=cm)
ax.contourf(xh3, yh3, np.ma.masked_array(wsp3, dbz3<20), levels=levs, vmin=0, vmax=30, cmap=cm)
ax.contourf(xh4, yh4, np.ma.masked_array(wsp4, dbz4<20), levels=levs, vmin=0, vmax=30, cmap=cm)
ax.contourf(xh5, yh5, np.ma.masked_array(wsp5, dbz5<20), levels=levs, vmin=0, vmax=30, cmap=cm)
ax.contourf(xh6, yh6, np.ma.masked_array(wsp6, dbz6<20), levels=levs, vmin=0, vmax=30, cmap=cm)
ax.contourf(xh7, yh7, np.ma.masked_array(wsp7, dbz7<20), levels=levs, vmin=0, vmax=30, cmap=cm)
ax.contourf(xh8, yh8, np.ma.masked_array(wsp8, dbz8<20), levels=levs, vmin=0, vmax=30, cmap=cm)

cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_ticks(np.linspace(0,30,7))
cb.set_label('Wind speed (m/s)', fontsize=10)

for i in range(8):
    is_rij = crit[f"{i+1}"]['is_rij']
    is_mv = crit[f"{i+1}"]['is_mv']
    is_db = crit[f"{i+1}"]['is_db']
    is_mv_rij = crit[f"{i+1}"]['is_mv_rij']
    is_mv_db = crit[f"{i+1}"]['is_mv_db']
    
    rij_mask = (is_rij<1) | (is_mv_rij)
    mv_mask = (is_mv<1) | (is_mv_rij) | (is_mv_db)
    db_mask = (is_db<1) | (is_mv_db)
    mv_rij_mask = (is_mv_rij<1)
    mv_db_mask = (is_mv_db<1)
    
    x = xh + 100 + i*x_added
    y = yh + 100 + i*y_added
    
    plot_cfill(x, y, np.ma.masked_array(is_rij, rij_mask), 'w', ax, datalims=[0,1], cmap='PiYG', cbar=False, alpha=0.6)
    plot_cfill(x, y, np.ma.masked_array(is_mv, mv_mask), 'w', ax, datalims=[0,1], cmap='bwr', cbar=False, alpha=0.6)
    plot_cfill(x, y, np.ma.masked_array(is_db, db_mask), 'w', ax, datalims=[0,1], cmap='Bu10', cbar=False, alpha=0.6)
    plot_cfill(x, y, np.ma.masked_array(is_mv_rij, mv_rij_mask), 'w', ax, datalims=[0,1], cmap='vanimo_r', cbar=False, alpha=0.7)
    plot_cfill(x, y, np.ma.masked_array(is_mv_db, mv_db_mask), 'w', ax, datalims=[0,1], cmap='managua_r', cbar=False, alpha=0.7)
    
    ax.contour(x, y, is_rij, levels=[0.1], colors='green', linewidths=1)

ax.contour(xh1, yh1, wsp1, levels=[25], colors='k', linewidths=[1])
ax.contour(xh2, yh2, wsp2, levels=[25], colors='k', linewidths=[1])
ax.contour(xh3, yh3, wsp3, levels=[25], colors='k', linewidths=[1])
ax.contour(xh4, yh4, wsp4, levels=[25], colors='k', linewidths=[1])
ax.contour(xh5, yh5, wsp5, levels=[25], colors='k', linewidths=[1])
ax.contour(xh6, yh6, wsp6, levels=[25], colors='k', linewidths=[1])
ax.contour(xh7, yh7, wsp7, levels=[25], colors='k', linewidths=[1])
ax.contour(xh8, yh8, wsp8, levels=[25], colors='k', linewidths=[1])




# ax.contour(xh1, yh1, shs1, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh2, yh2, shs2, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh3, yh3, shs3, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh4, yh4, shs4, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh5, yh5, shs5, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh6, yh6, shs6, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh7, yh7, shs7, levels=levs, colors=cols, linewidths=lws)
# ax.contour(xh8, yh8, shs8, levels=levs, colors=cols, linewidths=lws)
ax.set_xlim(xl)
ax.set_ylim(yl)
ax.set_xlabel('Translated x (km)', fontsize=10)
ax.set_ylabel('y (km)', fontsize=10)
ax.set_title(f"{sounding_str} - 80-m wind swaths + wind mechanisms (0-8 h)", fontsize=10)
# l1, = ax.plot([-2,-1], [-2,-1], 'gray', linewidth=0.75)
# l2, = ax.plot([-2,-1], [-2,-1], 'k', linewidth=1)
# # ax.legend(handles=[l1,l2], labels=['300 m2/s2','500 m2/s2'], loc='lower right', fontsize=9)
# ax.legend(handles=[l2], labels=['SVR'], loc='lower right', fontsize=9)

l1, = ax.plot([-2,-1], [-2,-1], 'green', linewidth=1)
l2, = ax.plot([-2,-1], [-2,-1], 'r', linewidth=1)
l3, = ax.plot([-2,-1], [-2,-1], 'dodgerblue', linewidth=1)
l4, = ax.plot([-2,-1], [-2,-1], 'violet', linewidth=1)
l5, = ax.plot([-2,-1], [-2,-1], 'gold', linewidth=1)
l6, = ax.plot([-2,-1], [-2,-1], 'k', linewidth=1)
ax.legend(handles=[l1,l2,l3,l4,l5,l6], labels=['RIJ','MV','DB','MV+RIJ','MV+DB','SVR'], loc='lower right', fontsize=7)

ax.text(xt[0], yt[0], '1 h', fontsize=9, fontweight='bold')
ax.text(xt[1], yt[1], '2 h', fontsize=9, fontweight='bold')
ax.text(xt[2], yt[2], '3 h', fontsize=9, fontweight='bold')
ax.text(xt[3], yt[3], '4 h', fontsize=9, fontweight='bold')
ax.text(xt[4], yt[4], '5 h', fontsize=9, fontweight='bold')
ax.text(xt[5], yt[5], '6 h', fontsize=9, fontweight='bold')
ax.text(xt[6], yt[6], '7 h', fontsize=9, fontweight='bold')
ax.text(xt[7], yt[7], '8 h', fontsize=9, fontweight='bold')

if figsave:
    plt.savefig(fp+'wind_mechanism_swath_ERA5-1_test7.png', dpi=300)






